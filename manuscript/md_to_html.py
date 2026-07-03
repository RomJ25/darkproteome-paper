"""
Standalone Markdown -> print-ready HTML for the manuscript (Python stdlib only; no installs).
Purpose: produce manuscript.html -- open in a browser and print to PDF for bioRxiv.
Deliberately conservative: handles headers, bold (**), italic (*...* only — NOT _ , so subscripts
like T_N / D_N / theta_N survive), inline `code`, tables, ordered/unordered lists, --- rules,
and paragraphs. Unicode math (alpha, rho, theta, superscripts) passes through as-is.
Run:  python3 manuscript/md_to_html.py   ->  manuscript/manuscript.html
"""
import re, html as _html
from pathlib import Path

SRC = Path("manuscript/manuscript.md")
OUT = Path("manuscript/manuscript.html")

def inline(text: str) -> str:
    # escape HTML special chars first (math uses < > &), then add our own tags
    t = _html.escape(text, quote=False)
    t = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t)      # bold
    t = re.sub(r"`([^`]+?)`", r"<code>\1</code>", t)              # inline code
    t = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"<em>\1</em>", t)  # *italic* (single-star only)
    return t

def main():
    lines = SRC.read_text(encoding="utf-8").split("\n")
    doc_title = next((l[2:].strip() for l in lines if l.startswith("# ")), "Manuscript")
    out, i, n = [], 0, len(lines)
    para = []
    def flush_para():
        if para:
            out.append("<p>" + inline(" ".join(para).strip()) + "</p>")
            para.clear()
    while i < n:
        line = lines[i]
        s = line.strip()
        if not s:
            flush_para(); i += 1; continue
        # header
        m = re.match(r"^(#{1,6})\s+(.*)$", s)
        if m:
            flush_para(); lvl = len(m.group(1))
            out.append(f"<h{lvl}>" + inline(m.group(2)) + f"</h{lvl}>"); i += 1; continue
        # horizontal rule
        if re.match(r"^(-{3,}|\*{3,}|_{3,})$", s):
            flush_para(); out.append("<hr/>"); i += 1; continue
        # table: consecutive lines starting with '|'
        if s.startswith("|"):
            flush_para(); tbl = []
            while i < n and lines[i].strip().startswith("|"):
                tbl.append(lines[i].strip()); i += 1
            def cells(row): return [c.strip() for c in row.strip().strip("|").split("|")]
            out.append("<table>")
            if tbl:
                out.append("<thead><tr>" + "".join(f"<th>{inline(c)}</th>" for c in cells(tbl[0])) + "</tr></thead>")
                body = tbl[2:] if len(tbl) > 1 and set(tbl[1].replace("|","").replace("-","").replace(":","").strip()) <= set() else tbl[1:]
                # robust: treat 2nd row as separator if it's all dashes/colons
                if len(tbl) > 1 and re.match(r"^\|?[\s:\-|]+\|?$", tbl[1]):
                    body = tbl[2:]
                else:
                    body = tbl[1:]
                out.append("<tbody>")
                for row in body:
                    out.append("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in cells(row)) + "</tr>")
                out.append("</tbody>")
            out.append("</table>"); continue
        # lists
        if re.match(r"^[-*+]\s+", s) or re.match(r"^\d+\.\s+", s):
            flush_para()
            ordered = bool(re.match(r"^\d+\.\s+", s))
            tag = "ol" if ordered else "ul"
            out.append(f"<{tag}>")
            while i < n:
                ls = lines[i].strip()
                mm = re.match(r"^(?:[-*+]|\d+\.)\s+(.*)$", ls)
                if not mm: break
                out.append("<li>" + inline(mm.group(1)) + "</li>"); i += 1
            out.append(f"</{tag}>"); continue
        # paragraph text
        para.append(s); i += 1
    flush_para()

    css = """
    body{font-family:Georgia,'Times New Roman',serif;max-width:820px;margin:40px auto;padding:0 24px;
         line-height:1.5;color:#111;font-size:11.5pt}
    h1{font-size:18pt;line-height:1.25;margin:0 0 6px} h2{font-size:13.5pt;margin:22px 0 6px;border-bottom:1px solid #ddd;padding-bottom:3px}
    h3{font-size:11.8pt;margin:16px 0 4px;font-style:italic;font-weight:600}
    p{margin:7px 0;text-align:justify} code{font-family:'SF Mono',Menlo,Consolas,monospace;font-size:10pt;background:#f4f4f4;padding:0 2px}
    table{border-collapse:collapse;width:100%;margin:12px 0;font-size:10pt}
    th,td{border:1px solid #bbb;padding:5px 8px;text-align:left;vertical-align:top} th{background:#f0f0f0}
    hr{border:none;border-top:1px solid #ccc;margin:18px 0} ol,ul{margin:7px 0 7px 22px} li{margin:3px 0}
    em{font-style:italic} strong{font-weight:700}
    @media print{body{margin:0;max-width:none;font-size:10.5pt} h2{page-break-after:avoid} table,figure{page-break-inside:avoid}}
    """
    html_doc = ("<!doctype html><html><head><meta charset='utf-8'>"
                f"<title>{_html.escape(doc_title)}</title>"
                f"<style>{css}</style></head><body>\n" + "\n".join(out) + "\n</body></html>")
    OUT.write_text(html_doc, encoding="utf-8")
    print(f"wrote {OUT}  ({len(html_doc):,} bytes, {len(out)} blocks)")
    print("Open it in a browser and Print -> Save as PDF for the bioRxiv manuscript file.")

if __name__ == "__main__":
    main()
