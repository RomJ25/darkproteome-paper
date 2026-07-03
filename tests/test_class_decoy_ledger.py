#!/usr/bin/env python3
"""Regression test for class_decoy_ledger.py (stdlib only; no pytest needed).

Runs the CLI on a synthetic mokapot/Percolator-style PSM table
(tests/sample_mokapot_psms.tsv) and checks the per-class ledger against the
hand-computed expectation. This exercises the --psms / q-value acceptance path —
the primary (mokapot) input mode, complementing the pepXML path demonstrated in
examples/. Run:  python3 tests/test_class_decoy_ledger.py
"""
import json, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TOOL = os.path.join(ROOT, "src", "darkproteome", "class_decoy_ledger.py")
FIXTURE = os.path.join(HERE, "sample_mokapot_psms.tsv")
MULT_FIXTURE = os.path.join(HERE, "sample_mokapot_psms_multiplicity.tsv")

# At alpha=0.01 (accept q<=0.01): rows c1-c3,d1,n1-n3,d2,v1 accepted; r1-r3 rejected.
EXPECTED_CLASS = {
    "canonical":    (3, 1),   # (T_class, D_class)
    "noncanonical": (3, 1),
    "variant":      (1, 0),
}
EXPECTED_GLOBAL = (7, 2)      # (global_T, global_D)


def main():
    out = os.path.join(tempfile.mkdtemp(), "ledger")
    cmd = [sys.executable, TOOL, "--psms", FIXTURE,
           "--id-col", "PSMId", "--decoy-col", "is_decoy",
           "--q-col", "mokapot_qvalue", "--accession-col", "Proteins",
           "--alpha", "0.01", "--out", out]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("FAIL — tool exited %d\n%s" % (r.returncode, r.stderr)); sys.exit(1)

    data = json.load(open(out + ".json"))
    led = {row["class"]: row for row in data["ledger"]}
    fails = []

    for cls, (T, D) in EXPECTED_CLASS.items():
        got = led.get(cls)
        if not got or (got["T_class"], got["D_class"]) != (T, D):
            fails.append("  class %s: expected T=%d D=%d, got %s" % (cls, T, D, got))
        elif abs(got["class_fdr_hat"] - (D + 1) / T) > 1e-4:
            fails.append("  class %s: (D+1)/T mismatch: %s" % (cls, got["class_fdr_hat"]))

    gT, gD = EXPECTED_GLOBAL
    if (data["meta"]["global_T"], data["meta"]["global_D"]) != (gT, gD):
        fails.append("  global: expected T=%d D=%d, got T=%s D=%s"
                     % (gT, gD, data["meta"]["global_T"], data["meta"]["global_D"]))
    if data["meta"]["fdr_basis"] != "reported q-values":
        fails.append("  fdr_basis: expected 'reported q-values', got %r" % data["meta"]["fdr_basis"])

    # --emit-diagfdr: the per-PSM diagFDR contract (id, is_decoy, q, pep, run, score)
    dfdr = out + ".diagfdr.tsv"
    r2 = subprocess.run(cmd + ["--emit-diagfdr", dfdr, "--run", "testrun"],
                        capture_output=True, text=True)
    if r2.returncode != 0:
        fails.append("  --emit-diagfdr exited %d: %s" % (r2.returncode, r2.stderr))
    else:
        dl = [ln.rstrip("\n").split("\t") for ln in open(dfdr)]
        if dl[0] != ["id", "is_decoy", "q", "pep", "run", "score"]:
            fails.append("  diagFDR header wrong: %s" % dl[0])
        if len(dl) - 1 != 12:
            fails.append("  diagFDR rows: expected 12, got %d" % (len(dl) - 1))
        ndec = sum(1 for row in dl[1:] if row[1] == "1")
        if ndec != 3:
            fails.append("  diagFDR is_decoy=1 count: expected 3, got %d" % ndec)

    if fails:
        print("FAIL — class_decoy_ledger:\n" + "\n".join(fails)); sys.exit(1)
    print("PASS — mokapot/--psms path: canonical 3T/1D, noncanonical 3T/1D, variant 1T/0D, "
          "global 7T/2D, (D+1)/T correct, q-value basis; --emit-diagfdr contract (6 cols, 12 PSMs, 3 decoys).")

    test_multiplicity()
    test_validation()


def test_validation():
    """Regression coverage for 4 previously-silent failure modes: each used to either
    silently produce wrong/empty output or crash with an unhelpful error; now each exits
    non-zero with a specific, actionable message."""
    fails = []

    def expect_error(name, cmd, needle):
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode == 0:
            fails.append(f"  {name}: expected a non-zero exit, tool ran to completion instead")
        elif needle not in r.stderr:
            fails.append(f"  {name}: expected {needle!r} in stderr, got: {r.stderr[-300:]}")

    out = os.path.join(tempfile.mkdtemp(), "ledger")
    base = [sys.executable, TOOL, "--psms", FIXTURE, "--id-col", "PSMId",
            "--decoy-col", "is_decoy", "--q-col", "mokapot_qvalue", "--alpha", "0.01"]

    # 1. --psms without --accession-col: used to silently classify every PSM "noncanonical".
    expect_error("missing --accession-col", base + ["--out", out], "--accession-col is required")

    # 2. --unit unique-peptide without --peptide-col: used to silently collapse every PSM
    #    of a (class, decoy-status) pair into one fake peptide.
    expect_error("missing --peptide-col for --unit unique-peptide",
                 base + ["--accession-col", "Proteins", "--unit", "unique-peptide", "--out", out],
                 "--peptide-col is required")

    # 3. A typo'd --q-col: used to silently produce an empty ledger (every PSM filtered out).
    bad_qcol = [sys.executable, TOOL, "--psms", FIXTURE, "--id-col", "PSMId",
                "--decoy-col", "is_decoy", "--q-col", "not_a_real_column",
                "--accession-col", "Proteins", "--alpha", "0.01", "--out", out]
    expect_error("typo'd --q-col", bad_qcol, "matched no values")

    # 4. A wrong --pepxml-score-name: used to silently parse 0 PSMs from a valid pepXML.
    pepxml = os.path.join(ROOT, "data", "external", "pxd055609_pepxml")
    if os.path.isdir(pepxml):
        files = [f for f in os.listdir(pepxml) if f.endswith(".pepXML")]
        if files:
            bad_score = [sys.executable, TOOL, "--pepxml", os.path.join(pepxml, files[0]),
                         "--pepxml-score-name", "not_a_real_score", "--alpha", "0.03", "--out", out]
            expect_error("wrong --pepxml-score-name", bad_score, "0 PSMs parsed")

    if fails:
        print("FAIL — class_decoy_ledger validation:\n" + "\n".join(fails)); sys.exit(1)
    print("PASS — 4 previously-silent failure modes now exit with a specific error "
          "(missing --accession-col/--peptide-col, typo'd --q-col, wrong --pepxml-score-name).")


# MULT_FIXTURE: canonical has CANONK (3 PSMs, target), OTHERCANONK (1 PSM, target),
# REVCANONK (1 PSM, decoy); noncanonical has CRYPTICK (2 PSMs, target), CRYPTICK2
# (1 PSM, target), REVCRYPTICK + REVCRYPTICK2 (1 PSM each, decoy). All q<=0.01 -> all accepted.
EXPECTED_UNIQUE_PEPTIDE = {
    "canonical":    (2, 1),   # {CANONK, OTHERCANONK} vs {REVCANONK}
    "noncanonical": (2, 2),   # {CRYPTICK, CRYPTICK2} vs {REVCRYPTICK, REVCRYPTICK2}
}
EXPECTED_STRATA = {
    "canonical":    {"n=1": (1, 1), "n=2": (0, 0), "n>=3": (1, 0)},
    "noncanonical": {"n=1": (1, 2), "n=2": (1, 0), "n>=3": (0, 0)},
}


def test_multiplicity():
    fails = []

    # --unit unique-peptide: dedupe by (class, peptide, decoy-status) before counting.
    out1 = os.path.join(tempfile.mkdtemp(), "ledger")
    cmd1 = [sys.executable, TOOL, "--psms", MULT_FIXTURE,
            "--id-col", "PSMId", "--decoy-col", "is_decoy", "--q-col", "mokapot_qvalue",
            "--accession-col", "Proteins", "--peptide-col", "peptide",
            "--alpha", "0.01", "--unit", "unique-peptide", "--out", out1]
    r1 = subprocess.run(cmd1, capture_output=True, text=True)
    if r1.returncode != 0:
        print("FAIL — --unit unique-peptide exited %d\n%s" % (r1.returncode, r1.stderr)); sys.exit(1)
    led = {row["class"]: row for row in json.load(open(out1 + ".json"))["ledger"]}
    for cls, (T, D) in EXPECTED_UNIQUE_PEPTIDE.items():
        got = led.get(cls)
        if not got or (got["T_class"], got["D_class"]) != (T, D):
            fails.append("  unit=unique-peptide, class %s: expected T=%d D=%d, got %s" % (cls, T, D, got))

    # --stratify-multiplicity: per-class PSM-replication-depth buckets, counted at PSM level.
    out2 = os.path.join(tempfile.mkdtemp(), "ledger")
    cmd2 = [sys.executable, TOOL, "--psms", MULT_FIXTURE,
            "--id-col", "PSMId", "--decoy-col", "is_decoy", "--q-col", "mokapot_qvalue",
            "--accession-col", "Proteins", "--peptide-col", "peptide",
            "--alpha", "0.01", "--stratify-multiplicity", "--out", out2]
    r2 = subprocess.run(cmd2, capture_output=True, text=True)
    if r2.returncode != 0:
        print("FAIL — --stratify-multiplicity exited %d\n%s" % (r2.returncode, r2.stderr)); sys.exit(1)
    strata = json.load(open(out2 + ".json"))["multiplicity_strata"]
    for cls, buckets in EXPECTED_STRATA.items():
        for b, (T, D) in buckets.items():
            got = strata.get(cls, {}).get(b, {"T": 0, "D": 0})
            if (got["T"], got["D"]) != (T, D):
                fails.append("  stratify, class %s bucket %s: expected T=%d D=%d, got %s" % (cls, b, T, D, got))

    if fails:
        print("FAIL — class_decoy_ledger multiplicity features:\n" + "\n".join(fails)); sys.exit(1)
    print("PASS — --unit unique-peptide dedupes correctly (canonical 2T/1D, noncanonical 2T/2D); "
          "--stratify-multiplicity buckets correctly (canonical n=1 1T/1D n>=3 1T/0D; "
          "noncanonical n=1 1T/2D n=2 1T/0D).")


if __name__ == "__main__":
    main()
