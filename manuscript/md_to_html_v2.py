"""
Standalone Markdown -> print-ready HTML for manuscript_v2 + supplement_v2 (Python stdlib only; no installs).
Adapted from md_to_html.py (which targets the withdrawn manuscript.md) -- same conservative parser,
retargeted at the two live v2 source files. Does not modify md_to_html.py or its output.
Run:  python3 manuscript/md_to_html_v2.py   ->  manuscript/manuscript_v2.html, manuscript/supplement_v2.html
"""
import base64, re, html as _html
from pathlib import Path

def inline(text: str) -> str:
    t = _html.escape(text, quote=False)
    # Code spans are pulled out to placeholders FIRST and restored last, so a bare `*` inside a code
    # span (e.g. `` `RPL*` ``) can never be mistaken by the later bold/italic passes for the start or
    # end of an emphasis run -- that mismatch is exactly what produced malformed nested tags before.
    spans = []
    def stash(m):
        spans.append(m.group(1))
        return f"\x00{len(spans) - 1}\x00"
    t = re.sub(r"`([^`]+?)`", stash, t)
    t = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t)
    t = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"<em>\1</em>", t)
    t = re.sub(r"\x00(\d+)\x00", lambda m: f"<code>{spans[int(m.group(1))]}</code>", t)
    return t

_IMG_EXT = {".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg", ".gif": "image/gif"}

def embed_image(rel_path: str, base_dir: Path) -> str:
    """Embed an image as a base64 data: URI so the built HTML is self-contained (no relative-path
    dependency on figures_v2/ sitting alongside the file -- matters once this is emailed or opened
    from elsewhere)."""
    p = (base_dir / rel_path).resolve()
    mime = _IMG_EXT.get(p.suffix.lower(), "application/octet-stream")
    data = base64.b64encode(p.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{data}"

def parse_blocks(lines: list, base_dir: Path) -> list:
    out_blocks, i, n = [], 0, len(lines)
    para = []
    def flush_para():
        if para:
            out_blocks.append("<p>" + inline(" ".join(para).strip()) + "</p>")
            para.clear()
    while i < n:
        line = lines[i]
        s = line.strip()
        if not s:
            flush_para(); i += 1; continue
        if s.startswith(">"):
            flush_para(); bq = []
            while i < n:
                ls = lines[i].strip()
                if not ls.startswith(">"):
                    break
                inner = ls[1:]
                if inner.startswith(" "):
                    inner = inner[1:]
                bq.append(inner); i += 1
            out_blocks.append("<blockquote>" + "".join(parse_blocks(bq, base_dir)) + "</blockquote>")
            continue
        m = re.match(r"^!\[([^\]]*)\]\(([^)]+)\)$", s)
        if m:
            flush_para()
            alt, src = m.group(1), m.group(2)
            uri = embed_image(src, base_dir)
            # Bold "Figure N." label, roman explanation -- matches the Table N. caption style
            # (a full-italic caption reads as one dense block; a bold label lets the eye anchor).
            cm = re.match(r"^((?:Figure|Supplementary Figure)\s+\S+\.)\s(.*)$", alt, re.S)
            cap_html = f"<strong>{inline(cm.group(1))}</strong> {inline(cm.group(2))}" if cm else inline(alt)
            cap = f"<figcaption>{cap_html}</figcaption>" if alt else ""
            out_blocks.append(f'<figure><img src="{uri}" alt="{_html.escape(alt, quote=True)}"/>{cap}</figure>')
            i += 1; continue
        m = re.match(r"^(#{1,6})\s+(.*)$", s)
        if m:
            flush_para(); lvl = len(m.group(1))
            out_blocks.append(f"<h{lvl}>" + inline(m.group(2)) + f"</h{lvl}>"); i += 1; continue
        if re.match(r"^(-{3,}|\*{3,}|_{3,})$", s):
            flush_para(); out_blocks.append("<hr/>"); i += 1; continue
        if s.startswith("|"):
            flush_para(); tbl = []
            while i < n and lines[i].strip().startswith("|"):
                tbl.append(lines[i].strip()); i += 1
            def cells(row): return [c.strip() for c in row.strip().strip("|").split("|")]
            out_blocks.append("<table>")
            if tbl:
                out_blocks.append("<thead><tr>" + "".join(f"<th>{inline(c)}</th>" for c in cells(tbl[0])) + "</tr></thead>")
                if len(tbl) > 1 and re.match(r"^\|?[\s:\-|]+\|?$", tbl[1]):
                    body = tbl[2:]
                else:
                    body = tbl[1:]
                out_blocks.append("<tbody>")
                for row in body:
                    out_blocks.append("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in cells(row)) + "</tr>")
                out_blocks.append("</tbody>")
            out_blocks.append("</table>"); continue
        if re.match(r"^[-*+]\s+", s) or re.match(r"^\d+\.\s+", s):
            flush_para()
            ordered = bool(re.match(r"^\d+\.\s+", s))
            tag = "ol" if ordered else "ul"
            items = []
            while i < n:
                ls = lines[i].strip()
                if not ls:
                    break
                mm = re.match(r"^(?:[-*+]|\d+\.)\s+(.*)$", ls)
                if mm:
                    items.append([mm.group(1)])
                elif items:
                    items[-1].append(ls)
                else:
                    break
                i += 1
            out_blocks.append(f"<{tag}>")
            for parts in items:
                out_blocks.append("<li>" + inline(" ".join(parts)) + "</li>")
            out_blocks.append(f"</{tag}>"); continue
        para.append(s); i += 1
    flush_para()
    return out_blocks

def linkify_dois(body_html: str) -> str:
    """Turn 'doi:X' into a clickable https://doi.org/X link, scoped to the References section only
    -- a generic inline linkifier would also catch stray 'doi:' text elsewhere (e.g. this docstring's
    own examples) and Chrome's print-to-pdf turns a real <a href> into a clickable PDF annotation."""
    start = body_html.find("<h2>References</h2>")
    if start == -1:
        return body_html
    next_h2 = body_html.find("<h2>", start + 1)
    end = next_h2 if next_h2 != -1 else len(body_html)
    before, refs, after = body_html[:start], body_html[start:end], body_html[end:]
    refs = re.sub(r"doi:(\S+?)\.(?=\s|<|$)",
                  lambda m: f'doi:<a href="https://doi.org/{m.group(1)}">{m.group(1)}</a>.',
                  refs)
    return before + refs + after

CSS = """
    body{font-family:Georgia,'Times New Roman',serif;max-width:820px;margin:40px auto;padding:0 24px;
         line-height:1.5;color:#111;font-size:11.5pt}
    h1{font-size:18pt;line-height:1.25;margin:0 0 6px} h2{font-size:13.5pt;margin:22px 0 6px;border-bottom:1px solid #ddd;padding-bottom:3px}
    h3{font-size:11.8pt;margin:16px 0 4px;font-style:italic;font-weight:600}
    p{margin:7px 0;text-align:justify} code{font-family:'SF Mono',Menlo,Consolas,monospace;font-size:10pt;background:#f4f4f4;padding:0 2px}
    table{border-collapse:collapse;width:100%;margin:12px 0;font-size:10pt}
    th,td{border:1px solid #bbb;padding:5px 8px;text-align:left;vertical-align:top} th{background:#f0f0f0}
    hr{border:none;border-top:1px solid #ccc;margin:18px 0} ol,ul{margin:7px 0 7px 22px} li{margin:3px 0}
    em{font-style:italic} strong{font-weight:700}
    blockquote{margin:14px 0;padding:1px 0 1px 16px;border-left:3px solid #ccc}
    blockquote>*:first-child{margin-top:0} blockquote>*:last-child{margin-bottom:0}
    figure{margin:16px 0;text-align:center} figure img{max-width:100%;height:auto;border:1px solid #ddd}
    figcaption{font-size:10pt;color:#333;margin-top:6px;text-align:left}
    @media print{body{margin:0 auto;max-width:442pt;font-size:11pt} h1,h2{page-break-after:avoid} table,figure{page-break-inside:avoid}}
    """

def wrap_html(body_html: str, title: str, lang: str = "en-US") -> str:
    """Wrap already-rendered block HTML in the shared page skeleton (doctype/head/CSS/body) --
    shared by convert() (per-file .html builds) and build_pdf.py (which needs one combined
    document instead of two, so a merge step can't strip Chrome's own PDF tagging)."""
    return (f"<!doctype html><html lang='{_html.escape(lang, quote=True)}'>"
            f"<head><meta charset='utf-8'><title>{_html.escape(title)}</title>"
            f"<style>{CSS}</style></head><body>\n{body_html}\n</body></html>")

def convert(src: Path, out: Path):
    lines = src.read_text(encoding="utf-8").split("\n")
    doc_title = next((l[2:].strip() for l in lines if l.startswith("# ")), "Manuscript")
    out_blocks = parse_blocks(lines, src.parent)
    body_html = linkify_dois("\n".join(out_blocks))
    html_doc = wrap_html(body_html, doc_title)
    out.write_text(html_doc, encoding="utf-8")
    print(f"wrote {out}  ({len(html_doc):,} bytes, {len(out_blocks)} blocks)")

def main():
    convert(Path("manuscript/manuscript_v2.md"), Path("manuscript/manuscript_v2.html"))
    convert(Path("manuscript/supplement_v2.md"), Path("manuscript/supplement_v2.html"))

if __name__ == "__main__":
    main()
