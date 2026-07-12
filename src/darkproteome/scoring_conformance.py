"""Guard: ONE scorer, and no field may mean two things.

The 2026-07-12 external review found three defects that were all the same defect -- a value was
derived in one place and consumed in another that meant something different by it:

  * `reference_model.py` scored `present(reported_fdr)` -- the field is POPULATED -- and its
    survivorship curve called that stage "source-translation substantiated". An FDR of 1.0
    counted exactly like an FDR of 0.001.
  * `min_peptide_len` was filled by every ingester with the HLA LIGAND length, and read by
    `axes.score_source_orf` as the TRYPTIC peptide length supporting the source protein.
  * `tumor_specificity_basis="broad-normal-panel"` was assigned to a normal-ligandome search, a
    GTEx RNA threshold, and a cohort inclusion criterion alike, and all three strict-passed.

Fixing three instances does not stop the fourth. These checks do. Run them in CI, and before
quoting any number:

    python3 src/darkproteome/scoring_conformance.py

Exit 0 = conformant. Exit 1 = something is scoring that should not be, or a threshold has
drifted, or a retired column has come back.
"""
import ast
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import evidence_dimensions as ed  # noqa: E402
import consensus_bar  # noqa: E402
import schema  # noqa: E402

# Fields that carry EVIDENCE. Comparing one of these against a threshold is scoring, and scoring
# is axes.py's job alone.
EVIDENCE_FIELDS = {
    "reported_fdr", "psm_qvalue", "periodicity_pct", "n_unique_peptides",
    "ligand_len", "source_pep_len", "tumor_specificity_modality",
    "tumor_specificity_scope", "validation_level",
}

# Columns retired on 2026-07-12 because each meant two things. They must not come back.
RETIRED_COLUMNS = {"min_peptide_len", "tumor_specificity_basis"}

# May read evidence fields for a legitimate reason that is NOT scoring.
ALLOWED = {
    "evidence_dimensions.py":  "THE scorer",
    "axes.py":                 "retired stub; raises on import",
    "consensus_bar.py":        "the protein-existence bar; thresholds asserted equal below",
    "schema.py":               "declares the columns",
    "ingest_cohorts.py":       "WRITES the fields from source tables",
    "ingest_atlases.py":       "WRITES the fields from source tables",
    "scoring_conformance.py":  "this file",
    "reference_model.py":      "reads *_reported for REUSABILITY; delegates all verdicts",
    "probe_report.py":         "tabulates validation_level for description, does not score",
    "robustness.py":           "tabulates validation_level for description, does not score",
}

COMPARISONS = (ast.Lt, ast.LtE, ast.Gt, ast.GtE)


def find_scorers():
    """Any module outside ALLOWED that compares an evidence field to a threshold."""
    problems = []
    for fn in sorted(os.listdir(HERE)):
        if not fn.endswith(".py") or fn in ALLOWED:
            continue
        path = os.path.join(HERE, fn)
        with open(path, encoding="utf-8") as fh:
            src = fh.read()
        try:
            tree = ast.parse(src, filename=fn)
        except SyntaxError as e:
            problems.append(f"{fn}: does not parse ({e})")
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.Compare):
                continue
            if not any(isinstance(op, COMPARISONS) for op in node.ops):
                continue
            names = {n.value for n in ast.walk(node) if isinstance(n, ast.Constant)
                     and isinstance(n.value, str)}
            hit = names & EVIDENCE_FIELDS
            if hit:
                problems.append(
                    f"{fn}:{node.lineno}: compares evidence field(s) {sorted(hit)} to a "
                    f"threshold. That is SCORING -- it belongs in evidence_dimensions.py.")
    return problems


def find_retired_columns():
    problems = []
    for fn in sorted(os.listdir(HERE)):
        if not fn.endswith(".py") or fn == "scoring_conformance.py":
            continue
        with open(os.path.join(HERE, fn), encoding="utf-8") as fh:
            for i, line in enumerate(fh, 1):
                code = line.split("#", 1)[0]
                for col in RETIRED_COLUMNS:
                    if f'"{col}"' in code or f"'{col}'" in code or f"{col}=" in code:
                        problems.append(
                            f"{fn}:{i}: retired column {col!r} is back. It meant two different "
                            f"things to two different readers; that is why it was split.")
    return problems


def check_thresholds():
    """axes.py and consensus_bar.py must not drift apart on the protein-existence bar."""
    t = consensus_bar.THRESHOLDS
    pairs = [
        ("max_protein_fdr", t["max_protein_fdr"], ed.MAX_SOURCE_FDR),
        ("min_unique_peptides", t["min_unique_peptides"], ed.MIN_SOURCE_PEPTIDES),
        ("min_source_pep_len", t["min_source_pep_len"], ed.MIN_SOURCE_PEP_LEN),
        ("min_periodicity_pct", t["min_periodicity_pct"], ed.MIN_PERIODICITY),
    ]
    return [f"threshold {k!r} drifted: consensus_bar={a} vs axes={b}"
            for k, a, b in pairs if a != b]


def check_schema():
    problems = []
    for col in RETIRED_COLUMNS:
        if col in schema.COLUMNS:
            problems.append(f"schema.COLUMNS still declares the retired column {col!r}")
    for col in ("ligand_len", "source_pep_len", "tumor_specificity_modality",
                "tumor_specificity_scope", "psm_qvalue", "accepted_decoys_in_class"):
        if col not in schema.COLUMNS:
            problems.append(f"schema.COLUMNS is missing the split column {col!r}")
    return problems


# THE POSITIVE CONTROLS. One per independent route to a strict-pass -- not one per axis.
#
# `source_orf` can strict-pass two ways: Ribo-seq periodicity (Path A) or protein-level FDR +
# unique peptides + tryptic length (Path B). The first version of this check used only
# SYNTHETIC_STRONG, which carries periodicity -- so it passed via Path A and never touched Path B.
# We then set MIN_SOURCE_PEP_LEN to 999, making Path B impossible, and the guard still said PASS.
#
# That matters more than it looks: Path B is exactly the route that is vacuous on the real corpus
# (0/307,318 claims report an FDR). An untested control on the route that carries the finding is
# no control at all. So every route gets its own witness.
POSITIVE_CONTROLS = {
    "SYNTHETIC_STRONG": "source_orf via Path A (Ribo-seq periodicity >= 70%)",
    "SYNTHETIC_STRONG_PROTEOMICS": "source_orf via Path B (protein FDR + n_peptides + tryptic len)",
}


def check_bar_is_satisfiable():
    """Every dimension, by every route, must be adjudicable-and-supported by SOMETHING.

    Source translation is un-adjudicable on the live corpus -- not one of 307,318 claims reports
    the statistic it needs. The obvious objection is that we rigged an impossible bar. These
    controls answer it: hand-built claims that report everything the criterion asks for MUST come
    out adjudicable and supporting on every dimension. If one ever stops, the bar has become
    unsatisfiable and no zero we report means anything.

    That is the difference between "no claim clears the bar" (a finding) and "no claim CAN clear
    the bar" (an artifact). This check is what keeps us on the right side of it.
    """
    import csv
    sample = os.path.join(os.path.dirname(os.path.dirname(HERE)), "data", "sample_claims.csv")
    if not os.path.exists(sample):
        return [f"positive controls missing: {sample}"]
    with open(sample, newline="", encoding="utf-8") as fh:
        rows = {r["peptide_sequence"]: r for r in csv.DictReader(fh)}

    # EVERY dimension, including class_fdr_reconstructible. An earlier version excused that one
    # from the controls on the grounds that no real source reports D_N -- which is precisely the
    # reasoning that would hide a structural zero. Asserting its reachability instead immediately
    # exposed a real bug: it returned `adjudicable=True, outcome='not-adjudicable'` for a claim
    # that reported everything, because the outcome had been hardcoded. Do not exempt a dimension
    # because you are confident about it.
    testable = list(ed.DIMENSIONS)

    problems = []
    for name, route in POSITIVE_CONTROLS.items():
        ctrl = rows.get(name)
        if not ctrl:
            problems.append(f"positive control {name} ({route}) is gone from data/sample_claims.csv")
            continue
        dims = ed.score_all(ctrl)
        for k in testable:
            d = dims[k]
            if not d["adjudicable"]:
                missing = [r for r in ed.RUNGS[:-1] if not d[r]]
                problems.append(
                    f"THE BAR IS UNSATISFIABLE via [{route}]: dimension {k!r} is NOT ADJUDICABLE "
                    f"on {name} (missing: {missing}), a claim that reports everything the "
                    f"criterion asks for. A zero here would measure the scorer, not the claims.")
            elif d["outcome"] != ed.SUPPORTS:
                problems.append(
                    f"{name}: dimension {k!r} is adjudicable but outcome={d['outcome']!r}, "
                    f"expected {ed.SUPPORTS!r}.")
    return problems


def check_unassayed_is_never_failure():
    """THE 2026-07-13 INVARIANT. Absence of evidence must never be recorded as evidence of absence.

    The old scorer returned `fail` for `MS-presented` on the immunogenicity axis -- scoring a claim
    that NOBODY EVER ASSAYED as an empirical immunogenicity failure. 307,274 of 307,318 rows are
    MS-presented, so it reported adjudicability as 307,278 when the number of claims actually
    carrying a human T-cell result is 2. The external reviewer caught it; we had shipped it.

    This asserts the fix directly: a claim with an MS observation and no assay must come out
    NOT ADJUDICABLE on the T-cell dimension -- never `contradicts`, never a failure. Same for a
    claim with no normal-tissue evidence at all.
    """
    problems = []
    unassayed = {
        "peptide_sequence": "SYNTHETICX",
        "evidence_types": "immunopeptidomics",
        "validation_level": "MS-presented",
        "hla_allele": "HLA-A*02:01", "ligand_len": "9",
    }
    d = ed.human_tcell_assay(unassayed)
    if d["adjudicable"]:
        problems.append("UNASSAYED SCORED AS ADJUDICABLE: an MS-presented peptide with no T-cell "
                        "assay must not be adjudicable on human_tcell_assay.")
    if d["outcome"] == ed.CONTRADICTS:
        problems.append("UNASSAYED SCORED AS AN EMPIRICAL NEGATIVE: 'nobody ran the assay' is not "
                        "'the assay was negative'. This is the exact defect the reviewer caught.")
    if ed.human_tcell_state(unassayed) != "not-assayed":
        problems.append(f"human_tcell_state should be 'not-assayed', got "
                        f"{ed.human_tcell_state(unassayed)!r}")

    no_normal = dict(unassayed, tumor_specificity_modality="not reported",
                     tumor_specificity_scope="not reported")
    n = ed.normal_presentation(no_normal)
    if n["adjudicable"] or n["outcome"] == ed.CONTRADICTS:
        problems.append("NO NORMAL-TISSUE EVIDENCE SCORED AS AN EMPIRICAL NEGATIVE: a claim with "
                        "no reported normal-tissue check is un-adjudicable, not 'not specific'.")
    return problems


def check_ladder_is_monotone():
    """The reporting ladder must be NESTED: rung k counts claims clearing rungs 1..k.

    A ladder that is not nested is not a ladder, and drawing one as a figure actively misleads.
    Before the rungs were made cumulative, `human_tcell_assay` read
        asserted 44 -> claim_linked 4 -> quantitative 4 -> modality_appropriate 307,316
    JUMPING BACK UP, because "not a mouse assay" is trivially true of every claim never assayed at
    all. Rendered as a heatmap, that put a reassuring dark column exactly where the evidence is
    absent. Caught by LOOKING AT THE FIGURE, not by reading the code that made it.
    """
    import csv
    sample = os.path.join(os.path.dirname(os.path.dirname(HERE)), "data", "sample_claims.csv")
    with open(sample, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    m = ed.matrix(rows)
    problems = []
    for k in ed.DIMENSIONS:
        v = [m[k][r] for r in ed.RUNGS]
        for i in range(len(v) - 1):
            if v[i] < v[i + 1]:
                problems.append(
                    f"LADDER NOT MONOTONE on {k!r}: {ed.RUNGS[i]}={v[i]} < "
                    f"{ed.RUNGS[i+1]}={v[i+1]}. Rungs must be cumulative, or a figure drawn from "
                    f"them will show a rise where the evidence is absent.")
    return problems


def main():
    problems = (find_scorers() + find_retired_columns()
                + check_thresholds() + check_schema() + check_bar_is_satisfiable()
                + check_unassayed_is_never_failure() + check_ladder_is_monotone())
    if problems:
        print("SCORING CONFORMANCE: FAIL\n")
        for p in problems:
            print(f"  - {p}")
        print(f"\n{len(problems)} problem(s).")
        return 1
    print("SCORING CONFORMANCE: PASS")
    print("  - evidence_dimensions.py is the only module scoring evidence fields")
    print("  - no retired (double-meaning) column has come back")
    print("  - consensus_bar.py thresholds still agree with evidence_dimensions.py")
    print(f"  - schema declares all {len(schema.COLUMNS)} columns, including the split ones")
    print(f"  - POSITIVE CONTROLS ({len(POSITIVE_CONTROLS)}, one per route): each is adjudicable and")
    print("    SUPPORTS on every testable dimension -- the bar is satisfiable, so a zero on the")
    print("    real corpus is a finding about the claims, not an artifact of the scorer")
    print("  - UNASSAYED IS NEVER FAILURE: an MS-presented peptide with no T-cell assay comes out")
    print("    NOT ADJUDICABLE, never as an empirical negative (the 2026-07-13 reviewer catch)")
    print("  - LADDER IS MONOTONE: rungs are cumulative, so a figure drawn from them cannot show")
    print("    a rise where the evidence is absent")
    return 0


if __name__ == "__main__":
    sys.exit(main())
