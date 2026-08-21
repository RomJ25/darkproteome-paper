"""
Build manuscript/manuscript_v2_full.pdf from the .md sources -- regenerates the per-file HTML
(manuscript_v2.html / supplement_v2.html, unchanged), then builds ONE combined HTML document (same
content, joined with a page break) and prints THAT once via the local Chrome binary, and finally
uses pypdf only to add metadata and a navigable outline (bookmarks) to that single PDF. Replaces
the ad hoc process README.md flags as missing ("that invocation was run ad hoc; if the PDFs need
rebuilding, the command needs reconstructing").

Why one combined print instead of print-each-then-merge: checked directly -- headless Chrome's
`--print-to-pdf` already emits a **tagged** PDF (pdfinfo reports `Tagged: yes`) with no extra flag
needed, and every figure's `<img>` already carries real alt text (its full caption -- see
`md_to_html_v2.py`'s `embed_image`) that flows into that tag structure. But `pypdf.PdfWriter.append()`
does NOT preserve that structure tree across a multi-document merge (confirmed: a single-document
pypdf clone+rewrite keeps `Tagged: yes`; a two-document `.append()` merge does not, in this pypdf
version). Printing one already-combined document sidesteps the merge entirely, so the tagging and
alt text Chrome already gives us for free survive into the shipped file.

Deliberately NOT using Playwright: it isn't installed here, and installing it (`pip install
playwright && playwright install chromium`) would be a new dependency for something the existing
ad hoc process already did with the Chrome.app binary already on this machine. subprocess + the
`--print-to-pdf` CLI flag reproduces that, unchanged.

Open items, not solved here:
  - A short-title running footer (beyond the bare page number `add_page_numbers()` now stamps).
    Chrome's `--print-to-pdf` CLI flag still has no header/footer support -- that would need the
    DevTools Protocol's `Page.printToPDF(displayHeaderFooter=True, footerTemplate=...)`, which needs
    a WebSocket client (`websocket-client` -- not installed) to drive. Not worth a new dependency for
    a running title alone; revisit if that tooling is ever installed for another reason.
  - PDF linearization ("fast web view"): `qpdf` is not installed either. If it's ever added to the
    toolchain, route the final file through `qpdf --linearize in.pdf out.pdf`.
  - A full accessibility structure tree beyond what Chrome already produces automatically is a
    materially bigger effort and out of scope here -- what's added here is specifically the
    page-navigation/hyperlink layer (bookmarks, DOI links, /Lang), on top of tags Chrome already
    gives us and this script no longer throws away.

Run:  python3 manuscript/build_pdf.py   ->  manuscript/manuscript_v2_full.pdf
"""
import re
import subprocess
import sys
from pathlib import Path

import pypdf
from pypdf.generic import NameObject, TextStringObject

REPO = Path(__file__).resolve().parent.parent
MDIR = REPO / "manuscript"
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
PAGE_BREAK = '<hr style="margin:28px 0 22px;border:none;border-top:1.5px solid #999">'

# (bookmark title, search needle for locating the page -- an ASCII-safe prefix of the real heading
# so text-extraction quirks around em-dashes/curly quotes can't cause a false miss, level)
MANUSCRIPT_OUTLINE = [
    ("Title", None, 0),
    ("Abstract", "Abstract", 0),
    ("Introduction", "Introduction", 0),
    ("Results", "Results", 0),
    ("R1", "R1. Most of IEAtlas's cancer-catalogued sequences also occur in canonical proteins", 1),
    # reorder: the consequence result now runs directly after R1, ahead of the
    # library/detection material -- R2 and R3 swapped CONTENT, not position, to answer "so what?"
    # immediately after the headline rather than several pages later.
    ("R2", "R2. The consequence is observable inside the resource", 1),
    ("R3", "R3. Latent canonical ambiguity differs enormously between libraries", 1),
    ("R4", "R4. The additive statistical problem, and the remedy the field already demonstrated", 1),
    ("R5", "R5. Two distinct problems; one reporting remedy; and what it costs", 1),
    ("Discussion", "Discussion", 0),
    # "Methods" the bare word is unsafe as a search needle -- it occurs 4 times as ordinary prose
    # inside the Discussion section ("...its Methods did not describe...") before the real heading,
    # confirmed by direct inspection. "Reference" is unique in the whole document (its own
    # subsection heading, immediately after "## Methods") and anchors both bookmarks correctly.
    ("Methods", "Reference", 0),
    ("Reference (Methods)", "Reference", 1),
    ("The overlap measurement", "The overlap measurement", 1),
    ("The library measurement", "The library measurement", 1),
    ("The detection-bias test", "The detection-bias test", 1),
    ("The normal-tissue consequence", "The normal-tissue consequence", 1),
    ("Class strata", "Class strata", 1),
    ("Class-specific FDR identifiability", "Class-specific FDR identifiability", 1),
    ("ImmunoVerse recurrence pilot", "ImmunoVerse recurrence and raw-lineage pilot", 1),
    ("Reproducibility", "Reproducibility", 1),
    # "References" the bare word is also unsafe: the Reproducibility subsection's last sentence
    # reads "...every entry in the References list below was independently re-verified..." --
    # currently on the same physical page as the real heading (so it happened to still resolve
    # correctly), but that's coincidental pagination, not a robust match. Anchor on the References
    # section's own header note instead, confirmed unique in the document.
    ("References", "Alphabetical by first author surname", 0),
]

SUPPLEMENT_OUTLINE = [
    # The bare word "Supplement" is unsafe as a search needle -- it occurs as ordinary prose several
    # times before the real heading (e.g. "...in the Supplement." in R4 and in "What this paper does
    # not claim"), confirmed by direct inspection to land the bookmark one page early. Anchor on the
    # Supplement's own opening sentence instead, confirmed unique in the document.
    ("Supplement", "Companion to the main manuscript", 0),
    ("S1", "S1. Where it can be interrogated, the ambiguity is structured by homology", 1),
    ("S2", "S2. Set-identification of a class-specific FDR", 1),
    ("S3", "S3. Recurrence and scan-level provenance in a frozen ImmunoVerse pilot", 1),
]


def normalize(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def find_pages(reader: "pypdf.PdfReader", needles: list) -> list:
    """For each needle (in document order), find the first page at or after a monotonically
    advancing pointer whose extracted text contains it -- so a repeated word can't jump a later
    bookmark backward. Falls back to a from-the-top search (and a loud warning) rather than
    silently mis-placing a bookmark."""
    page_texts = [normalize(p.extract_text() or "") for p in reader.pages]
    results, ptr = [], 0
    for needle in needles:
        if needle is None:
            results.append(ptr)
            continue
        n = normalize(needle)
        found = next((i for i in range(ptr, len(page_texts)) if n in page_texts[i]), None)
        if found is None:
            found = next((i for i in range(len(page_texts)) if n in page_texts[i]), None)
        if found is None:
            print(f"  ! WARNING: could not locate heading for bookmark: {needle!r} -- "
                  f"placing at page {ptr}")
            found = ptr
        results.append(found)
        ptr = found
    return results


def add_page_numbers(writer: "pypdf.PdfWriter"):
    """Stamp a centered page number near the bottom of every page -- without the DevTools Protocol
    header/footer machinery the module docstring rules out (needs websocket-client, a new dependency).

    First attempt was a per-page `pypdf.PageObject.merge_page()` of a Chrome-printed numbered overlay:
    technically correct output, but it more than DOUBLED the file (1.5 MB -> 3.6 MB) even though
    per-page content-stream size, font-object count and image bytes were all confirmed unchanged --
    merge_page's own page-merge machinery inflates the file by some other mechanism not worth chasing
    further once a clean alternative existed. That alternative, used here: append one small raw content
    stream directly to each page's own /Contents, referencing a single shared **base-14 Helvetica**
    font (Type1, no embedding -- every PDF viewer ships it, so this adds no font data at all). Measured
    cost: +4.7 KB total across 32 pages, vs +2 MB for the merge_page route.

    The stamp is PREPENDED to each page's /Contents array, not appended. A PDF content stream begins
    execution with an identity CTM by spec; Chrome's own content for the page may leave the graphics
    state transformed by the time it ends (confirmed: appending after it placed the number in the
    wrong spot, near the top-left instead of bottom-center). Prepending guarantees our absolute
    coordinates are interpreted before anything has transformed them, and our own `q ... Q` wrapper
    keeps the page's real content stream unaffected by our (identity) transform either way.

    Helvetica's digits are fixed-width (0.556 em each in the standard AFM metrics), so centering a
    1- or 2-digit page number needs no font-metrics library -- the width is exact, not estimated.

    Two follow-up fixes found by re-checking the shipped PDF at the content-stream level, not just by
    rendering it: (1) the stamp is wrapped in `/Artifact BMC ... EMC` so accessibility tools don't read
    folios into the reading order as body content. (2) the stamp's raw bytes must end with whitespace.
    Content streams in a page's /Contents array are logically concatenated with NO separator inserted
    -- per spec, applications must treat it as if it were one stream, byte for byte. Chrome's own
    per-page content immediately starts with a bare numeric operand for its global scale transform,
    e.g. `.23999999 0 0 -.23999999 0 792 cm`. Without a trailing space/newline, our closing `Q` glues
    directly onto that leading `.` -- pdfminer.six's tokenizer (unlike Poppler's, which is more
    permissive) reads `Q.23999999` as one unrecognized keyword, silently drops it, and the `cm` that
    follows fires with a corrupted operand stack. Confirmed live: before this fix, pdfplumber reported
    ~97% of ALL body-text characters on every page as having wildly out-of-page coordinates (e.g. y
    down to -2321 on a 792pt-tall page) -- not just the stamp, the entire page, because the botched
    global transform propagates to everything after it. Poppler (pdftotext/pdftoppm) tolerates the
    glued token and renders fine, which is exactly why this survived visual inspection undetected."""
    font_dict = pypdf.generic.DictionaryObject({
        NameObject("/Type"): NameObject("/Font"),
        NameObject("/Subtype"): NameObject("/Type1"),
        NameObject("/BaseFont"): NameObject("/Helvetica"),
        NameObject("/Encoding"): NameObject("/WinAnsiEncoding"),
    })
    font_ref = writer._add_object(font_dict)

    FONT_SIZE = 9
    DIGIT_WIDTH = 0.556 * FONT_SIZE
    BOTTOM_MARGIN = 28

    for i, page in enumerate(writer.pages, start=1):
        page_w = float(page.mediabox.width)
        text = str(i)
        x = (page_w - DIGIT_WIDTH * len(text)) / 2
        stream = pypdf.generic.StreamObject()
        stream.set_data(
            f"/Artifact BMC q 0.35 0.35 0.35 rg BT /PgNum {FONT_SIZE} Tf {x:.2f} {BOTTOM_MARGIN} Td "
            f"({text}) Tj ET Q EMC\n".encode("latin-1")
        )
        stream_ref = writer._add_object(stream)

        resources = page.get("/Resources")
        if resources is None:
            resources = pypdf.generic.DictionaryObject()
            page[NameObject("/Resources")] = resources
        fonts = resources.get("/Font")
        if fonts is None:
            fonts = pypdf.generic.DictionaryObject()
            resources[NameObject("/Font")] = fonts
        fonts[NameObject("/PgNum")] = font_ref

        contents = page.get("/Contents")
        if isinstance(contents, pypdf.generic.ArrayObject):
            contents.insert(0, stream_ref)
        elif contents is not None:
            page[NameObject("/Contents")] = pypdf.generic.ArrayObject([stream_ref, contents])
        else:
            page[NameObject("/Contents")] = pypdf.generic.ArrayObject([stream_ref])


def print_to_pdf(html_path: Path, pdf_path: Path):
    subprocess.run(
        [CHROME, "--headless", "--disable-gpu", "--no-pdf-header-footer",
         f"--print-to-pdf={pdf_path}", f"file://{html_path.resolve()}"],
        check=True, capture_output=True,
    )
    print(f"  wrote {pdf_path.name}")


def doc_title(md_path: Path) -> str:
    for line in md_path.read_text(encoding="utf-8").split("\n"):
        if line.startswith("# "):
            return line[2:].strip()
    return "Manuscript"


def main():
    sys.path.insert(0, str(MDIR))
    import md_to_html_v2 as m2

    print("Regenerating standalone HTML from .md sources...")
    m2.convert(MDIR / "manuscript_v2.md", MDIR / "manuscript_v2.html")
    m2.convert(MDIR / "supplement_v2.md", MDIR / "supplement_v2.html")

    print("Building one combined document (avoids a PDF-level merge, so Chrome's own tagging"
          " and figure alt text survive into the shipped file)...")
    title = doc_title(MDIR / "manuscript_v2.md")
    m_lines = (MDIR / "manuscript_v2.md").read_text(encoding="utf-8").split("\n")
    s_lines = (MDIR / "supplement_v2.md").read_text(encoding="utf-8").split("\n")
    m_body = m2.linkify_dois("\n".join(m2.parse_blocks(m_lines, MDIR)))
    s_body = "\n".join(m2.parse_blocks(s_lines, MDIR))
    combined_html = m2.wrap_html(m_body + PAGE_BREAK + s_body, title, lang="en-US")
    combined_path = MDIR / "_combined_full.html"
    combined_path.write_text(combined_html, encoding="utf-8")

    print("Printing to PDF via headless Chrome...")
    out = MDIR / "manuscript_v2_full.pdf"
    print_to_pdf(combined_path, out)
    combined_path.unlink()

    print("Locating bookmark targets...")
    reader = pypdf.PdfReader(out)
    outline = MANUSCRIPT_OUTLINE + SUPPLEMENT_OUTLINE
    pages = find_pages(reader, [t[1] for t in outline])

    print("Adding metadata, /Lang, and bookmarks (single-document rewrite -- no merge, so"
          " Chrome's tagging is preserved)...")
    writer = pypdf.PdfWriter(clone_from=out)
    print("Stamping page numbers (raw content-stream append, base-14 font, no new dependency)...")
    add_page_numbers(writer)
    writer.add_metadata({
        "/Title": title,
        "/Author": "Rom Jan",
        "/Subject": "Canonical-sequence overlap, provenance retention and confidence propagation "
                    "in public non-canonical HLA-peptide atlases",
        "/Keywords": "immunopeptidomics, ncORF, HLA, tumor antigen, canonical overlap, provenance, "
                     "aggregation, false discovery rate",
    })
    writer.root_object[NameObject("/Lang")] = TextStringObject("en-US")
    parent = None
    for (label, _needle, level), page in zip(outline, pages):
        if level == 0:
            parent = writer.add_outline_item(label, page)
        else:
            writer.add_outline_item(label, page, parent=parent)

    with open(out, "wb") as f:
        writer.write(f)
    print(f"wrote {out}  ({len(reader.pages)} pages)")


if __name__ == "__main__":
    main()
