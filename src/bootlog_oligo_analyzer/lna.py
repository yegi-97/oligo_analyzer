"""Nearest-neighbor Tm model for DNA oligos with LNA bases (LNA strand vs. DNA target).

DNA steps, per-end initiation, symmetry and salt correction follow primer3's default
method (SantaLucia 1998 unified parameters, SantaLucia salt correction, divalent cations
converted to monovalent equivalents), so plain DNA gives the same Tm as primer3.calc_tm
for oligos up to 60 nt. Steps that touch an LNA use the full
LNA.DNA parameters of Owczarzy et al. (2011) Biochemistry 50:9352, Table 2 (1 M Na+);
initiation for terminal LNAs is taken as identical to DNA, as in that paper.
"""

import math

# 5'->3' top strand step: (dH kcal/mol, dS cal/mol/K), SantaLucia 1998
DNA_NN = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
INIT_AT = (2.3, 4.1)  # per terminal A/T base pair
INIT_GC = (0.1, -2.8)  # per terminal G/C base pair
SYMMETRY_DS = -1.4

# "+X" = LNA nucleotide; key is the 5'->3' LNA-strand step. Owczarzy et al. 2011, Table 2.
LNA_NN = {
    "+A+A": (-9.991, -27.175),
    "+A+C": (-11.389, -28.963),
    "+A+G": (-12.793, -31.607),
    "+A+T": (-14.703, -40.75),
    "+C+A": (-14.177, -35.498),
    "+C+C": (-15.399, -36.375),
    "+C+G": (-14.558, -35.239),
    "+C+T": (-15.737, -41.218),
    "+G+A": (-13.959, -35.097),
    "+G+C": (-16.109, -40.738),
    "+G+G": (-13.022, -29.673),
    "+G+T": (-17.361, -45.858),
    "+T+A": (-10.318, -26.108),
    "+T+C": (-9.166, -21.535),
    "+T+G": (-10.046, -22.591),
    "+T+T": (-10.419, -27.683),
    "+AA": (-7.193, -19.723),
    "+AC": (-7.269, -18.336),
    "+AG": (-7.536, -18.387),
    "+AT": (-4.918, -12.943),
    "+CA": (-7.451, -18.38),
    "+CC": (-5.904, -11.904),
    "+CG": (-9.815, -23.491),
    "+CT": (-7.092, -16.825),
    "+GA": (-5.038, -11.656),
    "+GC": (-10.16, -24.651),
    "+GG": (-10.844, -26.58),
    "+GT": (-8.612, -22.327),
    "+TA": (-7.246, -19.738),
    "+TC": (-6.307, -15.515),
    "+TG": (-10.04, -25.744),
    "+TT": (-6.372, -16.902),
    "A+A": (-6.908, -18.135),
    "A+C": (-5.51, -11.824),
    "A+G": (-9.0, -22.826),
    "A+T": (-5.384, -13.537),
    "C+A": (-7.142, -18.333),
    "C+C": (-5.937, -12.335),
    "C+G": (-10.876, -27.918),
    "C+T": (-9.471, -25.07),
    "G+A": (-7.756, -19.302),
    "G+C": (-10.725, -25.511),
    "G+G": (-8.943, -20.833),
    "G+T": (-9.035, -22.742),
    "T+A": (-5.609, -16.019),
    "T+C": (-7.591, -19.031),
    "T+G": (-6.335, -15.537),
    "T+T": (-5.574, -14.149),
}

R = 1.987  # cal/mol/K
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def tm_nn(seq: str, lna: frozenset[int], mv_conc: float, dv_conc: float,
          dntp_conc: float, dna_conc: float) -> float:
    """Tm (C) of a concrete A/C/G/T sequence; `lna` holds 0-based LNA positions."""
    dh = ds = 0.0
    for i in range(len(seq) - 1):
        if i in lna or i + 1 in lna:
            key = ("+" if i in lna else "") + seq[i] + ("+" if i + 1 in lna else "") + seq[i + 1]
            h, s = LNA_NN[key]
        else:
            h, s = DNA_NN[seq[i:i + 2]]
        dh, ds = dh + h, ds + s
    for end in (seq[0], seq[-1]):
        h, s = INIT_AT if end in "AT" else INIT_GC
        dh, ds = dh + h, ds + s
    # Self-complementary: symmetry penalty and full strand concentration (as in primer3).
    symmetric = seq == seq[::-1].translate(_COMPLEMENT)
    if symmetric:
        ds += SYMMETRY_DS

    if dv_conc > 0 and dv_conc > dntp_conc:
        mv_conc += 120 * math.sqrt(dv_conc - dntp_conc)
    ds += 0.368 * (len(seq) - 1) * math.log(mv_conc / 1000)
    conc = dna_conc / (1e9 if symmetric else 4e9)
    return dh * 1000 / (ds + R * math.log(conc)) - 273.15
