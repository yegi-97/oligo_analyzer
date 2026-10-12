"""Explanatory notes, method references and conditions shown by each tab's Info button."""

from bootlog_oligo_analyzer import core

DEGENERATE = (
    "Degenerate sequences (IUPAC codes) are expanded into all concrete variants. Above "
    f"{core.MAX_TM_VARIANTS:,} (Tm), {core.MAX_HAIRPIN_VARIANTS:,} (hairpins) or "
    f"{core.MAX_DIMER_VARIANTS:,} (dimers) variant combinations, a reproducible random sample of "
    "that size is used instead; the detail view says so."
)

NAMES = ("Names are labels only. Every analysis runs on every sequence and every pair in input "
         "order; '#n' is the input position.")

INPUT = (
    "Input: 'name<TAB>sequence', 'name: sequence' (also ; or ,), 'name sequence', a bare "
    "sequence, a name line followed by the sequence on the next line, or FASTA (read with "
    "Bio.SeqIO; text before the first '>' and ';' lines are ignored). The name ends at the "
    "first tab, ':', ';' or ','. All whitespace inside a sequence line is ignored (spaces, "
    "tabs, non-breaking and zero-width spaces), so codon-style groups such as 'GGT GGT TTC' "
    "work; trailing tabs or spaces after a name are fine. A sequence must be on one line "
    "(except in FASTA). '+N' marks an LNA base, 5'-/-3' marks are stripped, U is read as T."
)


def conditions(cond: core.Conditions, with_temp: bool) -> str:
    text = (f"Conditions used: Na+/K+ {cond.mv_conc:g} mM, Mg2+ {cond.dv_conc:g} mM, "
            f"dNTP {cond.dntp_conc:g} mM, oligo {cond.dna_conc:g} nM")
    return text + (f", structures at {cond.temp_c:g} C." if with_temp else ".")


def properties(cond: core.Conditions) -> str:
    return "\n\n".join([
        "PROPERTIES",
        "Columns: Len = length, LNA = number of LNA bases, Variants = number of concrete "
        "sequences for IUPAC codes. Tm and GC columns show one value, or for degenerate "
        "sequences the mean with the (min-max) range over all variants. Click a row for all "
        "values, the reverse complement and the predicted hairpin.",
        "Tm primer3: primer3 calc_tm, SantaLucia 1998 nearest-neighbor parameters, SantaLucia "
        "salt correction with Mg2+ converted to Na+ equivalents (von Ahsen 2001). LNA bases are "
        "counted as DNA.",
        "Tm Biopython: Bio.SeqUtils.MeltingTemp.Tm_NN with the same Na+, Mg2+, dNTP and oligo "
        "concentration (dnac1 = dnac2 = oligo/2, the same C/4 term as primer3), nearest-neighbor "
        "table DNA_NN3 (Allawi & SantaLucia 1997) and salt correction method 7 (Owczarzy et al. "
        "2008, Biochemistry 47:5336). With Biopython's default salt correction both models are "
        "identical, so the difference between the two columns is the effect of the salt "
        "model. LNA bases are counted as DNA.",
        "Tm LNA: the primer3 model with every nearest-neighbor step that touches an LNA replaced "
        "by the LNA/DNA parameters of Owczarzy et al. 2011, Biochemistry 50:9352, Table 2. "
        "Above 60 nt the LNA shift is added to primer3's long-sequence Tm. The DNA salt "
        "correction is reused. '-' = no LNA in the sequence.",
        "GC %: exact min, mean and max over all variants (no sampling needed).",
        "Reverse complement (detail view): Bio.Seq.reverse_complement, IUPAC-aware, plain DNA.",
        DEGENERATE, NAMES, INPUT, conditions(cond, with_temp=False),
    ])


def _structures_common(cond: core.Conditions, limit: float) -> list[str]:
    return [
        "primer3 cannot model LNA. Hairpins and dimers are predicted on the plain DNA "
        "sequence, with LNA bases treated as DNA. In the drawing LNA bases have a purple ring.",
        f"dG at {cond.temp_c:g} C. WARN marks dG <= {limit} kcal/mol. For degenerate "
        "sequences the most stable variant (lowest dG) is shown; the detail view names it.",
        f"primer3 limit: hairpins need <= {core.THAL_MAX_LEN} nt; for dimers at least one "
        f"sequence must be <= {core.THAL_MAX_LEN} nt.",
        "Drawing: each base is a labelled circle, lines along a strand are the backbone and "
        "dark rungs connect paired bases. The first sequence is drawn in blue, the second in "
        "orange (outline, backbone and its 5'/3' labels); the detail text uses the same colours.",
        DEGENERATE, NAMES, conditions(cond, with_temp=True),
    ]


def hairpins(cond: core.Conditions) -> str:
    return "\n\n".join([
        "HAIRPINS",
        "primer3 calc_hairpin (primer3 thermodynamic alignment). "
        "Individual sequences: each sequence on its own. Joined pairs: every ordered pair "
        "joined 5'->3' (A + B and B + A), e.g. to check a fusion or ligation product. "
        "Click a row for the drawing.",
        "Hairpin drawing: the stem is two parallel strands (5' strand on the left going up, "
        "3' strand on the right coming down), the loop bases sit on a semicircle on top, "
        "unpaired ends hang below. Circles are coloured by base (A green, T red, G yellow, "
        "C blue).",
        *_structures_common(cond, core.HAIRPIN_WARN_DG),
    ])


def dimers(cond: core.Conditions) -> str:
    return "\n\n".join([
        "DIMERS",
        "primer3 calc_homodimer (each sequence with itself) and calc_heterodimer (every pair "
        "of sequences). Click a row for the drawing.",
        "Dimer drawing: the first sequence (blue) runs 5'->3' on top, the second (orange) "
        "3'->5' below, paired bases aligned in columns and joined by rungs; unpaired bases are "
        "shifted outward. Each strand is labelled with its name at its 5' start and has its own "
        "5' and 3' marks, so the top strand starts on the left and the bottom strand on the "
        "right. A homodimer shows two copies of the same sequence (strand 1 and strand 2).",
        "Detail text: each strand line, including the four alignment lines, starts with its "
        "sequence name and is printed in the strand colour.",
        *_structures_common(cond, core.DIMER_WARN_DG),
    ])
