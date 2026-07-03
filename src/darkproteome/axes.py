"""Four independent evidence axes for a dark-proteome tumor-antigen claim.

Why four axes and not one bar: a single HLA-I antigen is ONE short 8-12mer ligand, so the
protein-existence consensus bar (≥2 unique tryptic peptides; see consensus_bar.py) is the
right standard for the SOURCE-ORF axis only. A claim is a STRICT SURVIVOR only if it
strict-passes all four:

  1. source_orf        — is the source ORF plausibly translated?
  2. hla_presentation  — is the peptide actually presented on HLA?
  3. tumor_specificity — is it absent from normal presentation/expression?
  4. immunogenicity    — do T cells respond?

Verdicts: 'strict-pass' | 'weak-pass' | 'fail' | 'unverifiable'. Scored ONLY as reported.

Caveats (honest, surfaced not buried):
- hla_presentation strict-pass assumes class-specific FDR (the corpus can't show it from
  reported fields); a raw-spectra re-score is the downstream confirmation.
- tumor_specificity strict-pass means "no PUBLIC normal evidence reported", never "safe".
- immunogenicity strict-pass = T-cell-validated; in-vivo (humanized mouse) is weak-pass.
"""

NOT_REPORTED_TOKENS = {"", "not reported", "not stated", "na", "n/a", "nr", "none"}
AXES = ["source_orf", "hla_presentation", "tumor_specificity", "immunogenicity"]


def _num(value):
    if value is None:
        return None
    s = str(value).strip().lower()
    if s in NOT_REPORTED_TOKENS:
        return None
    s = s.rstrip("%")
    try:
        return float(s)
    except ValueError:
        return None


def _txt(value):
    return (value or "").strip().lower()


def score_source_orf(row):
    ev = _txt(row.get("evidence_types"))
    periodicity = _num(row.get("periodicity_pct"))
    fdr = _num(row.get("reported_fdr"))
    npep = _num(row.get("n_unique_peptides"))
    minlen = _num(row.get("min_peptide_len"))
    riboseq = "ribo" in ev
    # Path A: Ribo-seq translation with frame periodicity.
    if riboseq and periodicity is not None and periodicity >= 70:
        return "strict-pass"
    # Path B: protein-level proteomics evidence (the Prensner bar). NB: `min_peptide_len`
    # is overloaded — here it stands in for the protein-existence tryptic-peptide length
    # (>=9aa), while hla_presentation reads the SAME field as the HLA-ligand length (8-12).
    # This path requires fdr AND npep AND minlen all reported together, which fires on 0 rows
    # of the live corpus (atlases omit npep+fdr) — it is the bar for a hypothetical
    # fully-reported claim, not a scorer that currently bites.
    if fdr is not None and npep is not None and minlen is not None:
        if fdr <= 0.001 and npep >= 2 and minlen >= 9:
            return "strict-pass"
        return "fail"
    # Translation asserted (Ribo-seq) but not quantified.
    if riboseq:
        return "weak-pass"
    return "unverifiable"


def score_hla_presentation(row):
    ev = _txt(row.get("evidence_types"))
    minlen = _num(row.get("min_peptide_len"))
    hla = _txt(row.get("hla_allele"))
    eluted = "immunopeptidom" in ev or "eluted" in ev or "hla-ms" in ev
    proteomics = "proteomic" in ev or ("ms" in ev and not eluted)
    predicted_only = "predict" in ev and not (eluted or proteomics)
    if predicted_only:
        return "fail"  # predicted binder only, no observed presentation
    if eluted:
        good_len = minlen is not None and 8 <= minlen <= 12
        has_allele = hla not in NOT_REPORTED_TOKENS
        return "strict-pass" if (good_len and has_allele) else "weak-pass"
    if proteomics:
        return "weak-pass"  # MS for the protein, but not eluted-ligand evidence
    return "unverifiable"


def score_tumor_specificity(row):
    basis = _txt(row.get("tumor_specificity_basis"))
    if basis == "broad-normal-panel":
        # strongest reported basis (claim checked a broad normal panel) -> strict-pass.
        # Still only "no public normal evidence was reported", NOT "proven safe".
        return "strict-pass"
    if basis == "matched-normal-only":
        return "weak-pass"
    return "unverifiable"


def score_immunogenicity(row):
    level = _txt(row.get("validation_level"))
    if level == "t-cell-validated":
        return "strict-pass"
    if level == "in-vivo":
        return "weak-pass"     # humanized-mouse < autologous human T cells
    if level == "validated_negative":
        return "fail"          # assayed and NOT reactive (a benchmark negative)
    if level in ("ms-presented", "nominated"):
        return "fail"          # presented/nominated, no immunogenicity evidence
    return "unverifiable"


def score_all(row):
    return {
        "source_orf": score_source_orf(row),
        "hla_presentation": score_hla_presentation(row),
        "tumor_specificity": score_tumor_specificity(row),
        "immunogenicity": score_immunogenicity(row),
    }


def is_strict_survivor(verdicts):
    return all(v == "strict-pass" for v in verdicts.values())
