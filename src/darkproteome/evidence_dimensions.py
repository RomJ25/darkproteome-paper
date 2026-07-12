"""Evidence DIMENSIONS: the reporting-and-adjudicability matrix.

The auditor scores a published claim on six evidence dimensions, and never as a single pass/fail
survivor count. A joint pass/fail across dimensions the reported record cannot even decide is a
measurement of the scorer, not of the claims.


THE REPORTING LADDER (per claim i, per dimension k)
---------------------------------------------------
A study can do excellent science and still fail the claim-linked reusability level. That level is
this audit's subject, so the ladder is the unit of analysis. The rungs are CUMULATIVE -- rung k
counts claims that clear rungs 1..k:

  asserted             the experiment/analysis is named for this claim
  claim_linked         an individual result travels WITH the claim (not just a study-level method)
  quantitative         that result is a value a prespecified criterion can be applied to
  modality_appropriate the measurement answers the endpoint actually being claimed
  adjudicable          all of the above -> the criterion can be applied independently

  outcome  (meaningful only when adjudicable)
      supports | contradicts | indirect | not-adjudicable

Two rules govern every dimension, and violating either manufactures results:

  A ZERO IS NOT A FINDING UNTIL THE RECORD COULD HAVE SAID OTHERWISE.
      A structural zero (no claim reports the field the criterion needs, so none CAN clear it) and
      an empirical zero (claims report it, the criterion applies, none clear it) look identical in
      a survivor count and mean opposite things. `adjudicable` separates them.

  ABSENCE OF EVIDENCE IS NEVER EVIDENCE OF ABSENCE.
      A peptide observed by MS with no T-cell assay is NOT a negative immunogenicity result: it is
      unassayed. A ligand length outside the HLA-I window is not a contradiction: no source here
      states the MHC class, and 13-25mers are ordinary HLA-II ligands. Both are `not-adjudicable`.

`adjudicable` is deliberately NOT named `A`. Source-attribution resolution -- whether the reported
record leaves only the nominated source compatible, i.e. O_i(R,E) = {u_i*} -- is a different
object: a record can carry excellent claim-linked translation evidence for the nominated ORF while
a canonical source remains perfectly compatible. It is computed from sequence-compatibility,
elsewhere, and never from this module.


WHY SIX DIMENSIONS AND NOT FOUR
-------------------------------
Presentation, allele restriction and search-error reconstructibility are separate questions, and
conjoining them hides coverage collapses. A single "HLA presentation" test requiring an eluted
ligand AND an assigned allele reports one verdict for two very different facts: elution is
asserted almost everywhere, while the allele the peptide was restricted to is reported for a tiny
minority. The conjunction reads as good coverage; the split shows the truth.
"""
DIMENSIONS = [
    "source_translation",         # is the nominated ORF translated?
    "hla_elution",                # was the peptide observed on HLA by MS?
    "allele_restriction",         # WHICH allele presented it?
    "normal_presentation",        # is it absent from NORMAL HLA presentation?
    "human_tcell_assay",          # do human T cells respond?
    "class_fdr_reconstructible",  # can the class-specific identification error be recomputed?
]

NOT_REPORTED_TOKENS = {"", "not reported", "not stated", "na", "n/a", "nr", "none",
                       "insufficient-info"}

# Prespecified criteria. Named once; `scoring_conformance.py` asserts consensus_bar.py agrees.
MAX_SOURCE_FDR = 0.001
MIN_SOURCE_PEPTIDES = 2
MIN_SOURCE_PEP_LEN = 9
MIN_PERIODICITY = 70.0
LIGAND_LEN_RANGE = (8, 12)

MODALITY_LIGANDOME_BROAD = "normal-ligandome-broad"
MODALITY_LIGANDOME_MATCHED = "normal-ligandome-matched"
MODALITY_RNA_BROAD = "normal-rna-broad"
MODALITY_RNA_MATCHED = "normal-rna-matched"
LIGANDOME_MODALITIES = (MODALITY_LIGANDOME_BROAD, MODALITY_LIGANDOME_MATCHED)
RNA_MODALITIES = (MODALITY_RNA_BROAD, MODALITY_RNA_MATCHED)
ALL_MODALITIES = LIGANDOME_MODALITIES + RNA_MODALITIES

SCOPE_PER_CLAIM = "per-claim-reported"
SCOPE_INCLUSION = "cohort-inclusion-criterion"

SUPPORTS, CONTRADICTS, INDIRECT, NOT_ADJ = "supports", "contradicts", "indirect", "not-adjudicable"


def _num(v):
    if v is None:
        return None
    s = str(v).strip().lower()
    if s in NOT_REPORTED_TOKENS:
        return None
    try:
        return float(s.rstrip("%"))
    except ValueError:
        return None


def _txt(v):
    return (v or "").strip().lower()


def _rep(v):
    return _txt(v) not in NOT_REPORTED_TOKENS


def _rung(asserted=False, claim_linked=False, quantitative=False,
          modality_appropriate=False, outcome=NOT_ADJ):
    adjudicable = bool(asserted and claim_linked and quantitative and modality_appropriate)
    return {"asserted": bool(asserted), "claim_linked": bool(claim_linked),
            "quantitative": bool(quantitative),
            "modality_appropriate": bool(modality_appropriate),
            "adjudicable": adjudicable,
            "outcome": outcome if adjudicable else NOT_ADJ}


# ---------------------------------------------------------------------------------------
# 1. source_translation
# ---------------------------------------------------------------------------------------
# Three levels that must be kept apart:
#   (a) translation ASSERTED, and QUANTIFIED AT THE STUDY LEVEL. The audited sources DO report
#       thresholds -- a RibORF score cutoff, an average read periodicity, a PSM-level FDR. It is
#       simply WRONG to say "the field does not report translation statistics"; it reports them for
#       the STUDY.
#   (b) CLAIM-LINKED quantitative support -- the individual ncORF's own score, the individual
#       peptide's own q-value. This is what a reader would need to re-adjudicate ONE claim, and it
#       is what is not published.
#   (c) translation biologically absent -- NOT ASSESSED, and not assessable from a reported record.
# Collapsing (a) into (b) mischaracterises what the authors did, and that error is fatal to an audit.
# The ladder encodes the distinction: `asserted` passes, `claim_linked` fails for the statistic, so
# `quantitative` and `adjudicable` are 0.

def source_translation(row):
    ev = _txt(row.get("evidence_types"))
    riboseq = "ribo" in ev
    ms = "immunopeptidom" in ev or "proteomic" in ev or "ms" in ev
    asserted = riboseq or ms
    # A binary RibORF call IS a claim-linked result -- it just is not a statistic.
    claim_linked = asserted
    periodicity = _num(row.get("periodicity_pct"))
    fdr = _num(row.get("reported_fdr"))
    npep = _num(row.get("n_unique_peptides"))
    slen = _num(row.get("source_pep_len"))
    quant_ribo = riboseq and periodicity is not None
    quant_prot = fdr is not None and npep is not None and slen is not None
    quantitative = quant_ribo or quant_prot
    outcome = NOT_ADJ
    if quant_ribo:
        outcome = SUPPORTS if periodicity >= MIN_PERIODICITY else CONTRADICTS
    elif quant_prot:
        outcome = (SUPPORTS if (fdr <= MAX_SOURCE_FDR and npep >= MIN_SOURCE_PEPTIDES
                                and slen >= MIN_SOURCE_PEP_LEN) else CONTRADICTS)
    return _rung(asserted, claim_linked, quantitative, asserted, outcome)


# ---------------------------------------------------------------------------------------
# 2. hla_elution   |   3. allele_restriction
# ---------------------------------------------------------------------------------------
# Kept apart deliberately. Conjoining them lets near-universal ELUTION mask the near-total absence
# of ALLELE assignment: the largest cancer-epitope atlas in the corpus ships four columns
# (Sequence, Length, ORF_ID, tissue) and no allele at all.

def hla_elution(row):
    ev = _txt(row.get("evidence_types"))
    eluted = "immunopeptidom" in ev or "eluted" in ev or "hla-ms" in ev
    predicted_only = "predict" in ev and not eluted
    ln = _num(row.get("ligand_len"))
    seq = _txt(row.get("peptide_sequence"))
    claim_linked = eluted and _rep(seq)
    quantitative = ln is not None
    lo, hi = LIGAND_LEN_RANGE
    outcome = NOT_ADJ
    if quantitative:
        # A length outside the HLA-I window is NOT a contradiction: none of these sources states
        # the MHC class, and 13-25mers are ordinary HLA-II ligands. Calling them `contradicts`
        # would be the same sin as scoring an unassayed claim as a failure -- inventing a negative
        # the record never gave. They are INDIRECT: consistent with presentation, not on the
        # criterion we prespecified.
        outcome = SUPPORTS if lo <= ln <= hi else INDIRECT
    if predicted_only:
        # a predicted binder is not an observation of presentation, at any modality
        return _rung(True, claim_linked, quantitative, False, NOT_ADJ)
    return _rung(eluted, claim_linked, quantitative, eluted, outcome)


def allele_restriction(row):
    ev = _txt(row.get("evidence_types"))
    eluted = "immunopeptidom" in ev or "eluted" in ev or "hla-ms" in ev
    has_allele = _rep(row.get("hla_allele"))
    # `quantitative` reads as "a specific, structured value the criterion applies to" -- here, a
    # named allele rather than a free-text gesture at HLA typing.
    return _rung(asserted=eluted, claim_linked=has_allele, quantitative=has_allele,
                 modality_appropriate=eluted, outcome=SUPPORTS if has_allele else NOT_ADJ)


# ---------------------------------------------------------------------------------------
# 4. normal_presentation
# ---------------------------------------------------------------------------------------
# Zero strict-passes here is NOT "these claims are biologically non-specific". It is "no claim
# carries modality-appropriate, claim-linked evidence of absence from the NORMAL LIGANDOME".
# Three states, kept distinct:
#   modality-appropriate + claim-linked  -> adjudicable
#   indirect / study-level only          -> RNA evidence, or a ligandome used only as an
#                                           inclusion criterion (cannot be re-derived per claim)
#   explicit normal DETECTION            -> the only true empirical negative for tumour absence.
#                                           Not an as-reported field anywhere in this corpus; it
#                                           is an AUTHOR-COMPUTED overlap and is reported
#                                           separately (never mixed into as-reported columns).

RESULT_ABSENT = "absent-from-normal"      # the record says: not found in the queried normals
RESULT_DETECTED = "detected-in-normal"    # the record says: IT WAS FOUND. an empirical negative.


def normal_presentation(row):
    modality = _txt(row.get("tumor_specificity_modality"))
    scope = _txt(row.get("tumor_specificity_scope"))
    result = _txt(row.get("tumor_specificity_result"))
    asserted = modality in ALL_MODALITIES
    claim_linked = scope == SCOPE_PER_CLAIM and result in (RESULT_ABSENT, RESULT_DETECTED)
    ligandome = modality in LIGANDOME_MODALITIES
    # The criterion here is CATEGORICAL -- "does the record state this peptide was, or was not,
    # found in the queried normal tissue?" -- so a structured per-claim verdict IS the value the
    # criterion applies to. It must NOT be hardcoded False on the grounds that the underlying
    # FPKM/TPM is unavailable: that would make the dimension un-adjudicable by construction, i.e.
    # a structural zero introduced by the scorer itself.
    quantitative = claim_linked
    outcome = NOT_ADJ
    if result == RESULT_ABSENT:
        outcome = SUPPORTS
    elif result == RESULT_DETECTED:
        outcome = CONTRADICTS       # explicit normal detection: the one true empirical negative
    return _rung(asserted, claim_linked, quantitative, ligandome, outcome)


def normal_presentation_state(row):
    """`0 adjudicable` means NO claim carries modality-appropriate, claim-linked evidence of
    absence from the normal ligandome — NOT that the claims are biologically non-specific."""
    modality = _txt(row.get("tumor_specificity_modality"))
    scope = _txt(row.get("tumor_specificity_scope"))
    result = _txt(row.get("tumor_specificity_result"))
    if result == RESULT_DETECTED:
        return "explicit-normal-detection"      # empirical negative for tumour absence
    if modality not in ALL_MODALITIES:
        return "not-reported"
    if modality in LIGANDOME_MODALITIES and scope == SCOPE_PER_CLAIM:
        return "modality-appropriate-claim-linked"
    return "indirect-or-study-level"


# ---------------------------------------------------------------------------------------
# 5. human_tcell_assay
# ---------------------------------------------------------------------------------------
# `MS-presented` must NEVER score as a failure. A peptide observed by MS with no T-cell assay is
# not a negative immunogenicity result -- it was never assayed. Nearly the whole corpus is
# MS-presented, so treating that as a decided failure would report adjudicability for essentially
# every claim while the number carrying an actual human T-cell result is two.
#
# Note the fourth state, which no binary could express: a cohort may ASSAY a peptide and publish
# the per-peptide result only inside a figure. Assay asserted; result not claim-linked.

def human_tcell_assay(row):
    level = _txt(row.get("validation_level"))
    ev = _txt(row.get("evidence_types"))
    assay_run = level in ("t-cell-validated", "validated_negative", "in-vivo") or "t-cell" in ev
    if level == "t-cell-validated":
        return _rung(True, True, True, True, SUPPORTS)      # human assay, positive
    if level == "validated_negative":
        return _rung(True, True, True, True, CONTRADICTS)   # human assay, negative
    if level == "in-vivo":
        # humanized mouse: a real assay, but not the human endpoint being claimed
        return _rung(True, True, True, False, NOT_ADJ)
    # assayed, but the per-claim result lives only in a figure -> asserted, NOT claim-linked
    if assay_run:
        return _rung(True, False, False, True, NOT_ADJ)
    # ms-presented / nominated / not reported == NOT ASSAYED. Never a failure.
    return _rung(False, False, False, True, NOT_ADJ)


def human_tcell_state(row):
    level = _txt(row.get("validation_level"))
    ev = _txt(row.get("evidence_types"))
    if level == "t-cell-validated":
        return "human-assay-positive"
    if level == "validated_negative":
        return "human-assay-negative"
    if level == "in-vivo":
        return "nonhuman-assay-indirect"
    if "t-cell" in ev:
        return "assayed-result-not-claim-linked"
    return "not-assayed"


# ---------------------------------------------------------------------------------------
# 6. class_fdr_reconstructible
# ---------------------------------------------------------------------------------------
# Search-error reconstructibility is its own dimension, not a footnote to presentation.
#
# A class-specific FDR needs the per-class ACCEPTED DECOY count D_N. HCC's sheet S26 -- the only
# table in the entire corpus reporting per-PSM statistics -- carries a `target_decoy` column and
# 43 rows, EVERY ONE a target. The authors had the labels; the published rows are targets only.
# So a claim can be claim-linked AND quantitative here (a per-PSM q-value) and STILL not
# adjudicable, because D_N is absent. That is exactly the point, and this dimension is the only
# place in the matrix where those two things come apart.

def class_fdr_reconstructible(row):
    q = _num(row.get("psm_qvalue"))
    ev = _txt(row.get("evidence_types"))
    asserted = "immunopeptidom" in ev or "proteomic" in ev or "ms" in ev   # a TD search happened
    claim_linked = q is not None
    quantitative = q is not None
    # D_N -- the per-class accepted-decoy count -- is not reported anywhere in this corpus, so the
    # class-specific error rate cannot be reconstructed for ANY claim, however good its q-value.
    d_n_available = _rep(row.get("accepted_decoys_in_class"))
    # The outcome must NOT be hardcoded to NOT_ADJ: this dimension asks "is the class-specific
    # error rate reconstructible?", so a record that DOES report D_N supports it. Hardcoding the
    # outcome made a fully-reported claim come back `adjudicable=True, outcome='not-adjudicable'`
    # -- self-contradictory, and it would have silently mis-scored the first source ever to
    # publish D_N. Found by asserting reachability rather than assuming it.
    return _rung(asserted, claim_linked, quantitative, d_n_available,
                 SUPPORTS if d_n_available else NOT_ADJ)


_DIMS = {
    "source_translation": source_translation,
    "hla_elution": hla_elution,
    "allele_restriction": allele_restriction,
    "normal_presentation": normal_presentation,
    "human_tcell_assay": human_tcell_assay,
    "class_fdr_reconstructible": class_fdr_reconstructible,
}

RUNGS = ["asserted", "claim_linked", "quantitative", "modality_appropriate", "adjudicable"]


def score_all(row):
    return {k: _DIMS[k](row) for k in DIMENSIONS}


def matrix(rows):
    """The reporting-and-adjudicability matrix. Counts, per dimension, per rung.

    THE RUNGS ARE CUMULATIVE: rung k counts claims that clear rungs 1..k. A ladder that is not
    nested is not a ladder, and drawing one as a figure actively misleads.

    Counting each flag independently breaks this: "not a mouse assay" is trivially true of every
    claim never assayed at all, so an independently-counted `modality_appropriate` rung RISES above
    the `claim_linked` rung below it. Drawn as a heatmap that puts a reassuring dark column exactly
    where the evidence is absent. Cumulative counts can only decrease, so a figure built from them
    cannot show a rise where the record is silent.
    """
    out = {k: {r: 0 for r in RUNGS} for k in DIMENSIONS}
    for k in DIMENSIONS:
        out[k].update({"n": 0, SUPPORTS: 0, CONTRADICTS: 0})
    for row in rows:
        for k in DIMENSIONS:
            d = _DIMS[k](row)
            out[k]["n"] += 1
            alive = True
            for r in RUNGS:
                alive = alive and bool(d[r])
                if alive:
                    out[k][r] += 1
            if d["outcome"] in (SUPPORTS, CONTRADICTS):
                out[k][d["outcome"]] += 1
    return out


def format_matrix(m, title=""):
    w = 27
    lines = []
    if title:
        lines.append(title)
    lines.append(f"  {'evidence dimension':<{w}}" + "".join(f"{r.replace('_',' '):>15}" for r in RUNGS)
                 + f"{'supports':>11}")
    lines.append("  " + "-" * (w + 15 * len(RUNGS) + 11))
    for k in DIMENSIONS:
        d = m[k]
        cells = "".join(f"{d[r]:>15,}" for r in RUNGS)
        lines.append(f"  {k:<{w}}{cells}{d[SUPPORTS]:>11,}")
    lines.append(f"  {'(of n claims)':<{w}}{m[DIMENSIONS[0]]['n']:>15,}")
    return "\n".join(lines)
