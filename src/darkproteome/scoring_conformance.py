"""Conformance guard: one scorer, one meaning per field, and a bar that can actually be cleared.

    python3 src/darkproteome/scoring_conformance.py

Exit 0 = conformant. Exit 1 = something is scoring that should not be, a threshold has drifted, a
double-meaning column has returned, the evidence bar has become unsatisfiable, an unassayed claim
is being scored as a failure, or the reporting ladder is no longer nested.

Every defect this guards against is the same defect: A VALUE PRODUCED IN ONE PLACE AND CONSUMED
SOMEWHERE THAT MEANS SOMETHING DIFFERENT BY IT. Fixing instances does not stop the next one; these
checks do, so run them before quoting any number.

  ONE SCORER. `evidence_dimensions.py` alone compares an evidence field to a threshold. Anything
  else that does so is a second scorer, and second scorers drift. Enforced at the AST level.

  NO DOUBLE-MEANING COLUMNS. An HLA ligand length and the tryptic peptide length supporting a
  source protein are different measurements; a per-PSM q-value and a protein-level FDR are
  different units. Fields that once carried two meanings must not return.

  REPORTING-COMPLETENESS IS NOT EVIDENCE-STRENGTH. "The FDR field is populated" is not "the FDR
  clears the bar" -- an FDR of 1.0 populates the field exactly as well as an FDR of 0.001.

  THE BAR MUST BE SATISFIABLE. If a criterion cannot be cleared by anything, a zero measures the
  scorer rather than the literature. `data/sample_claims.csv` carries hand-built claims that
  report everything the criteria ask for; they must come out adjudicable and supporting on EVERY
  dimension. One witness per independent ROUTE to a pass -- a control that clears a dimension by
  one route never tests the other. And no dimension is exempt: exempting one because "no real
  source reports that field anyway" is exactly the reasoning that hides a structural zero.

  ABSENCE OF EVIDENCE IS NEVER EVIDENCE OF ABSENCE. An MS-presented peptide with no T-cell assay
  must come out NOT ADJUDICABLE -- never as an empirical negative. "Nobody ran the assay" is not
  "the assay was negative".

  THE LADDER MUST BE NESTED. Rung k counts claims clearing rungs 1..k. A ladder that is not
  cumulative can rise where the evidence is absent, and a figure drawn from it will say so.
"""
import ast
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import evidence_dimensions as ed  # noqa: E402
import consensus_bar  # noqa: E402
import schema  # noqa: E402

# Fields that carry EVIDENCE. Comparing one of these against a threshold is SCORING, and scoring
# belongs to evidence_dimensions.py alone.
EVIDENCE_FIELDS = {
    "reported_fdr", "psm_qvalue", "periodicity_pct", "n_unique_peptides",
    "ligand_len", "source_pep_len", "tumor_specificity_modality",
    "tumor_specificity_scope", "validation_level",
}

# Columns retired because each carried two different measurements. They must not come back.
RETIRED_COLUMNS = {"min_peptide_len", "tumor_specificity_basis"}

# May read evidence fields for a legitimate reason that is NOT scoring.
ALLOWED = {
    "evidence_dimensions.py":  "THE scorer",
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
    """evidence_dimensions.py and consensus_bar.py must not drift apart on the protein-existence bar."""
    t = consensus_bar.THRESHOLDS
    pairs = [
        ("max_protein_fdr", t["max_protein_fdr"], ed.MAX_SOURCE_FDR),
        ("min_unique_peptides", t["min_unique_peptides"], ed.MIN_SOURCE_PEPTIDES),
        ("min_source_pep_len", t["min_source_pep_len"], ed.MIN_SOURCE_PEP_LEN),
        ("min_periodicity_pct", t["min_periodicity_pct"], ed.MIN_PERIODICITY),
    ]
    return [f"threshold {k!r} drifted: consensus_bar={a} vs evidence_dimensions={b}"
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


# THE POSITIVE CONTROLS. One per independent ROUTE to a pass -- not one per dimension.
#
# Source translation can be satisfied two ways: Ribo-seq periodicity (Path A), or protein-level
# FDR + unique peptides + tryptic length (Path B). A single control carrying periodicity clears
# the dimension via Path A and never exercises Path B at all -- so Path B could be made impossible
# and the guard would still pass. Path B is the route that is un-adjudicable on the real corpus,
# i.e. the one that carries the finding. An untested control on that route is no control.
# Every route gets its own witness.
POSITIVE_CONTROLS = {
    "SYNTHETIC_STRONG": "source_orf via Path A (Ribo-seq periodicity >= 70%)",
    "SYNTHETIC_STRONG_PROTEOMICS": "source_orf via Path B (protein FDR + n_peptides + tryptic len)",
}


def check_bar_is_satisfiable():
    """Every dimension, by every route, must be adjudicable-and-supported by SOMETHING.

    Source translation is un-adjudicable on the audited corpus: no claim reports the statistic the
    criterion needs. The obvious objection is that the bar was rigged to be impossible. These
    controls answer it -- hand-built claims that report everything the criteria ask for MUST come
    out adjudicable and supporting on every dimension.

    That is the difference between "no claim clears the bar" (a finding about the literature) and
    "no claim CAN clear the bar" (an artifact of the scorer).
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
    """Absence of evidence must never be recorded as evidence of absence.

    Scoring `MS-presented` as a FAILURE on the immunogenicity dimension treats a claim that was
    never assayed as an empirical negative. Nearly the whole corpus is MS-presented, so that would
    report adjudicability for essentially every claim while the number carrying an actual human
    T-cell result is two.

    A claim with an MS observation and no assay must come out NOT ADJUDICABLE -- never
    `contradicts`, never a failure. Likewise a claim with no normal-tissue evidence at all.
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
                        "'the assay was negative'.")
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
    If the rungs are counted independently, `modality_appropriate` on the T-cell dimension RISES
    above the `claim_linked` rung below it -- because "not a mouse assay" is trivially true of
    every claim never assayed at all. Rendered as a heatmap that puts a reassuring dark column
    exactly where the evidence is absent.
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
    print("    NOT ADJUDICABLE, never as an empirical negative")
    print("  - LADDER IS MONOTONE: rungs are cumulative, so a figure drawn from them cannot show")
    print("    a rise where the evidence is absent")
    return 0


if __name__ == "__main__":
    sys.exit(main())
