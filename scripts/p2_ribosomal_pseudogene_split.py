"""Does P2's ribosomal-enrichment excess survive once ribosomal-NAMED PSEUDOGENES are separated
from TRUE ribosomal protein-coding genes?

WHY THIS SCRIPT EXISTS. An external review of manuscript_v2
pointed out that P2's `is_ribosomal(g)` test (abundance_bias.py) flags any gene symbol starting
with RPL/RPS/MRPL/MRPS -- and that this prefix rule cannot distinguish a TRUE ribosomal
protein-coding gene (e.g. RPS3A) from a ribosomal-NAMED PSEUDOGENE of one (e.g. RPS3AP12, the
manuscript's own worked example of pseudogene-driven canonical overlap, S1/supplement). Pseudogenes
are already independently established by this project (S1, pseudogene_parent_authoritative.py) to be
mechanically canonical-overlapping for a reason that has NOTHING to do with protein abundance or
detection bias -- they are retro-transposed near-copies of their parent's coding sequence. If
ribosomal-named pseudogenes are riding inside P2's "ribosomal" bucket, some or all of the reported
"excess" (the catalogue's ribosomal enrichment OVER the library's own composition baseline) could be
this already-known pseudogene mechanism wearing a different label, not new evidence of abundance-
driven overdetection. P2's existing design tests library-COMPOSITION vs. detection; it was never
built to separate true-ribosomal-gene from ribosomal-pseudogene within the "ribosomal" bucket itself.

This script re-derives P2 (catalogue-side and library-side ribosomal risk ratios) with the ribosomal
bucket split in two, using the AUTHORITATIVE NCBI gene_group pseudogene->parent registry already on
disk for S1 (data/external/pseudogene_parents/), not a symbol-string heuristic -- so a gene symbol
starting with RPL/RPS/MRPL/MRPS is called a "pseudogene" here only if NCBI's curated relation says so.

    python3 scripts/p2_ribosomal_pseudogene_split.py
"""
import csv
import json
import math
import os
import re
import sys
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXT = os.path.join(REPO, "data", "external")
IE = os.path.join(EXT, "atlases", "IEAtlas_Epitopes_In_Cancer_Tissues.txt")
NUORF = os.path.join(EXT, "nuorfdb", "PA_nuORFdb_v1.2_protein.fasta")
SCORED = os.path.join(REPO, "data", "claim_catalog_scored.csv")
csv.field_size_limit(10_000_000)
K = 9

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pseudogene_parent_authoritative as ppa  # noqa: E402 -- reuse the authoritative registry loader


def is_ribosomal(g):
    g = (g or "").upper()
    return g.startswith(("RPL", "RPS", "MRPL", "MRPS"))


def two_prop(k1, n1, k2, n2):
    p1, p2 = (k1 / n1 if n1 else float("nan")), (k2 / n2 if n2 else float("nan"))
    return p1, p2, (p1 / p2 if p2 else float("inf"))


def fasta(path):
    name, seq = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def main():
    for p in (IE, SCORED, NUORF, ppa.GENE_INFO, ppa.GENE_GROUP):
        if not os.path.exists(p):
            sys.exit(f"missing {p}")

    print("Loading NCBI's authoritative pseudogene->parent registry (same source as S1)...")
    sym2id, syn2id, id2sym, _ = ppa.load_gene_info()
    p2par = ppa.load_gene_group()
    pseudo_ids = set(p2par.keys())
    print(f"  {len(pseudo_ids):,} curated human pseudogene GeneIDs (NCBI gene_group)")

    def resolve_id(sym):
        s = (sym or "").upper()
        if s in sym2id:
            return sym2id[s]
        cands = syn2id.get(s)
        if cands and len(cands) == 1:
            return next(iter(cands))
        return None

    _cache = {}

    def category(sym):
        """-> 'pseudogene' | 'true_gene' | 'unresolved', for a ribosomal-prefixed symbol."""
        if sym in _cache:
            return _cache[sym]
        gid = resolve_id(sym)
        if gid is None:
            cat = "unresolved"
        elif gid in pseudo_ids:
            cat = "pseudogene"
        else:
            cat = "true_gene"
        _cache[sym] = cat
        return cat

    # ---------------------------------------------------------------- catalogue side (IEAtlas)
    selfmap = {r["peptide"]: int(r["canonical_self"])
               for r in csv.DictReader(open(SCORED, newline="", encoding="utf-8"))}

    tissues, gene = defaultdict(set), {}
    with open(IE, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for r in rd:
            if len(r) < 4 or not r[0]:
                continue
            p = r[0].strip().upper()
            if not p.isalpha():
                continue
            tissues[p].add((r[3] or "").strip())
            gene.setdefault(p, (r[2] or "").split("_")[0].upper())

    scored = [p for p in tissues if p in selfmap]
    ov = [p for p in scored if selfmap[p]]
    nov = [p for p in scored if not selfmap[p]]

    ribo_syms = sorted({gene[p] for p in scored if is_ribosomal(gene.get(p))})
    cats = {s: category(s) for s in ribo_syms}
    print(f"\nCatalogue: {len(ribo_syms)} distinct ribosomal-prefixed gene symbols "
          f"({sum(1 for c in cats.values() if c == 'pseudogene')} authoritative pseudogenes, "
          f"{sum(1 for c in cats.values() if c == 'true_gene')} true genes, "
          f"{sum(1 for c in cats.values() if c == 'unresolved')} unresolved)")
    print("  pseudogene symbols :", ", ".join(s for s in ribo_syms if cats[s] == "pseudogene") or "(none)")

    print("\n" + "=" * 92)
    print("P2, SPLIT: true ribosomal genes vs. ribosomal-NAMED pseudogenes (NCBI-authoritative)")
    print("=" * 92)
    print("  IN THE CATALOGUE:")
    cat_results = {}
    for label in ("true_gene", "pseudogene", "unresolved"):
        ro = sum(1 for p in ov if cats.get(gene.get(p)) == label)
        rn = sum(1 for p in nov if cats.get(gene.get(p)) == label)
        c1, c2, c_rr = two_prop(ro, len(ov), rn, len(nov))
        cat_results[label] = c_rr
        print(f"    [{label:>10}]  overlapping {ro:>6,}/{len(ov):<9,} = {100*c1:6.3f}%   "
              f"non-overlapping {rn:>6,}/{len(nov):<9,} = {100*c2:6.3f}%   RR = {c_rr:.2f}x")

    # ---------------------------------------------------------------- library side (nuORFdb)
    print("\n  IN THE LIBRARY (nuORFdb, distinct 9-mers) -- scanning...")
    SPROT = os.path.join(EXT, "swissprot_human.fasta")
    canon = set()
    for _, s in fasta(SPROT):
        for i in range(len(s) - K + 1):
            canon.add(s[i:i + K])

    kmers_by_cat = {"true_gene": set(), "pseudogene": set(), "unresolved": set(), "nonribo": set()}
    n_orf = 0
    orf_gene_cache = {}
    for h, s in fasta(NUORF):
        if len(s) < K:
            continue
        n_orf += 1
        m = re.search(r"GN=(\S+)", h)
        g = (m.group(1) if m else "").upper()
        if is_ribosomal(g):
            if g not in orf_gene_cache:
                orf_gene_cache[g] = category(g)
            lab = orf_gene_cache[g]
        else:
            lab = "nonribo"
        for i in range(len(s) - K + 1):
            kmers_by_cat[lab].add(s[i:i + K])
    print(f"  scanned {n_orf:,} nuORFdb ORFs")

    all_kmers = set().union(*kmers_by_cat.values())
    noncanon_kmers = all_kmers - canon
    canon_kmers = all_kmers & canon
    lib_results = {}
    for label in ("true_gene", "pseudogene", "unresolved"):
        rk = kmers_by_cat[label]
        ov_r = len(rk & canon_kmers)
        no_r = len(rk & noncanon_kmers)
        # Same convention as abundance_bias.py's P2: rate = (this subcategory's kmers) / (all
        # library kmers), evaluated separately within the canonical-overlapping and non-overlapping
        # partitions of the library.
        r1 = ov_r / len(canon_kmers) if canon_kmers else float("nan")
        r2 = no_r / len(noncanon_kmers) if noncanon_kmers else float("nan")
        rr = r1 / r2 if r2 else float("inf")
        lib_results[label] = rr
        print(f"    [{label:>10}]  {len(rk):>7,} distinct 9-mers   canonical-overlapping rate "
              f"{100*r1:6.3f}%   non-overlapping rate {100*r2:6.3f}%   RR = {rr:.2f}x")

    print("\n" + "=" * 92)
    print("EXCESS (catalogue RR / library RR), per sub-category -- this is the number that matters")
    print("=" * 92)
    for label in ("true_gene", "pseudogene", "unresolved"):
        c_rr, l_rr = cat_results[label], lib_results[label]
        excess = c_rr / l_rr if l_rr else float("inf")
        print(f"  [{label:>10}]  catalogue {c_rr:.2f}x / library {l_rr:.2f}x  =  excess {excess:.2f}x")

    excess = {label: round(cat_results[label] / lib_results[label], 2) if lib_results[label] else None
              for label in ("true_gene", "pseudogene", "unresolved")}
    out = {
        "ribo_symbols_total": len(ribo_syms),
        "ribo_symbols_pseudogene": [s for s in ribo_syms if cats[s] == "pseudogene"],
        "ribo_symbols_true_gene": [s for s in ribo_syms if cats[s] == "true_gene"],
        "ribo_symbols_unresolved": [s for s in ribo_syms if cats[s] == "unresolved"],
        "catalogue_rr": {k: round(v, 2) for k, v in cat_results.items()},
        "library_rr": {k: round(v, 2) for k, v in lib_results.items()},
        "excess": excess,
    }
    art = os.path.join(REPO, "data", "derived_p2_pseudogene_split.json")
    json.dump(out, open(art, "w"), indent=2)
    print(f"\nwrote {os.path.relpath(art, REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
