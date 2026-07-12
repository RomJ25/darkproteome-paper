"""The field's consensus evidence bar, as code.

Reference: Prensner et al., Mol Cell Proteomics 2023; HUPO-HPP audit (Wright et al. 2024).

We evaluate a claim ONLY against what it reports. We never recompute FDR from raw spectra
— that is a separate, downstream reanalysis. A claim that omits a required field is
`insufficient-info`, which is itself the most common — and most damning — outcome.
"""

try:
    from .schema import NOT_REPORTED
except ImportError:  # allow running as a top-level script (see audit.py)
    from schema import NOT_REPORTED

# Mirrors evidence_dimensions.MAX_SOURCE_FDR / MIN_SOURCE_PEPTIDES / MIN_SOURCE_PEP_LEN /
# MIN_PERIODICITY. scoring_conformance.py fails if these drift apart.
THRESHOLDS = {
    "max_protein_fdr": 0.001,     # <= 0.1%
    "min_unique_peptides": 2,
    "min_source_pep_len": 9,      # aa -- the TRYPTIC peptide supporting the source protein.
                                  # NOT the HLA ligand length (`ligand_len`): a different
                                  # measurement, and most HLA ligands are >= 9aa anyway, so
                                  # sharing one field would satisfy this rule by coincidence.
    "min_periodicity_pct": 70.0,
}

_RIBOSEQ_MARKERS = ("ribo", "riboseq", "ribo-seq", "periodicity")


def _num(value):
    """Parse a reported numeric field; return None for 'not reported'/blank/garbage."""
    if value is None:
        return None
    s = str(value).strip().lower()
    if s in ("", NOT_REPORTED, "na", "n/a", "nr", "none"):
        return None
    s = s.rstrip("%")
    try:
        return float(s)
    except ValueError:
        return None


def _evidence_uses_riboseq(evidence_types):
    s = (evidence_types or "").lower()
    return any(m in s for m in _RIBOSEQ_MARKERS)


def evaluate_consensus_bar(row):
    """Evaluate one claim against the bar, as reported.

    Returns {verdict, criteria} where verdict is 'yes' | 'no' | 'insufficient-info' and
    criteria maps each rule to 'pass' | 'fail' | 'unknown' | 'n/a'.
    """
    fdr = _num(row.get("reported_fdr"))
    n_pep = _num(row.get("n_unique_peptides"))
    min_len = _num(row.get("source_pep_len"))
    periodicity = _num(row.get("periodicity_pct"))
    needs_periodicity = _evidence_uses_riboseq(row.get("evidence_types"))

    criteria = {
        "protein_fdr<=0.1%": _check(fdr, lambda v: v <= THRESHOLDS["max_protein_fdr"]),
        ">=2_unique_peptides": _check(n_pep, lambda v: v >= THRESHOLDS["min_unique_peptides"]),
        "source_pep_len>=9": _check(min_len, lambda v: v >= THRESHOLDS["min_source_pep_len"]),
    }
    if needs_periodicity:
        criteria["riboseq_periodicity>=70%"] = _check(
            periodicity, lambda v: v >= THRESHOLDS["min_periodicity_pct"]
        )
    else:
        criteria["riboseq_periodicity>=70%"] = "n/a"

    required = [v for v in criteria.values() if v != "n/a"]
    if any(v == "unknown" for v in required):
        verdict = "insufficient-info"
    elif all(v == "pass" for v in required):
        verdict = "yes"
    else:
        verdict = "no"

    return {"verdict": verdict, "criteria": criteria}


def _check(value, predicate):
    if value is None:
        return "unknown"
    return "pass" if predicate(value) else "fail"
