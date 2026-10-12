"""Oligo calculations (Tm, GC, reverse complement, hairpins, dimers) based on primer3 and Biopython.

Names are labels only: every analysis runs on every sequence (and every pair) in input order,
and nothing is ever matched, paired or grouped by name.
"""

import io
import itertools
import math
import random
from collections.abc import Iterator
from dataclasses import dataclass, field

import primer3
from Bio import SeqIO
from Bio.Seq import reverse_complement
from Bio.SeqUtils import MeltingTemp

from bootlog_oligo_analyzer.lna import tm_nn

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

# Variants are enumerated exhaustively up to these limits, above them randomly sampled.
MAX_TM_VARIANTS = 50_000
MAX_HAIRPIN_VARIANTS = 100_000
MAX_DIMER_VARIANTS = 20_000
THAL_MAX_LEN = 60  # primer3 thermodynamic alignment limit

HAIRPIN_WARN_DG = -2.0  # kcal/mol
DIMER_WARN_DG = -6.0  # kcal/mol

# Second Tm model for comparison. With its defaults Tm_NN equals primer3 exactly (same
# nearest-neighbor set and salt correction), so it uses the Owczarzy 2008 Mg2+-aware salt
# correction (method 7) instead.
BIO_NN_TABLE = MeltingTemp.DNA_NN3  # Allawi & SantaLucia 1997 (Biopython default)
BIO_SALTCORR = 7  # Owczarzy et al. 2008


def ascii_only(text: str) -> str:
    return text.encode("ascii", "replace").decode("ascii")


@dataclass
class Conditions:
    mv_conc: float = 50.0  # monovalent cations, mM
    dv_conc: float = 1.5  # divalent cations, mM
    dntp_conc: float = 0.6  # mM
    dna_conc: float = 250.0  # oligo, nM
    temp_c: float = 37.0  # temperature for structure dG

    def tm_kwargs(self) -> dict:
        return dict(mv_conc=self.mv_conc, dv_conc=self.dv_conc,
                    dntp_conc=self.dntp_conc, dna_conc=self.dna_conc)

    def thermo_kwargs(self) -> dict:
        return dict(self.tm_kwargs(), temp_c=self.temp_c)


@dataclass(frozen=True)
class Oligo:
    number: int  # 1-based input position; identifies the sequence in reports
    name: str  # free-text label, never interpreted
    seq: str  # plain IUPAC sequence without LNA marks
    lna: frozenset[int] = frozenset()  # 0-based LNA positions

    @property
    def label(self) -> str:
        return f"#{self.number} {self.name}"

    @property
    def marked(self) -> str:
        """Sequence in +N notation."""
        return "".join(("+" if i in self.lna else "") + b for i, b in enumerate(self.seq))

    @property
    def n_variants(self) -> int:
        return math.prod(len(IUPAC[c]) for c in self.seq)


@dataclass
class MinMeanMax:
    min: float
    mean: float
    max: float


@dataclass
class BasicResult:
    oligo: Oligo
    tm_primer3: MinMeanMax
    tm_bio: MinMeanMax
    tm_lna: MinMeanMax | None  # None without LNA
    tm_checked: int
    tm_exhaustive: bool
    gc: MinMeanMax
    revcomp: str


@dataclass
class StructResult:
    oligos: tuple[Oligo, ...]  # 1 (hairpin, homodimer) or 2 (joined hairpin, heterodimer)
    kind: str  # "hairpin" or "dimer"
    found: bool = False
    dg: float = 0.0  # kcal/mol
    dh: float = 0.0  # kcal/mol
    ds: float = 0.0  # cal/mol/K
    tm: float = 0.0
    structure: list[str] = field(default_factory=list)  # primer3 ASCII lines, tab prefix removed
    variant: tuple[str, ...] = ()  # most stable concrete variant(s)
    n_total: int = 1
    n_checked: int = 1
    exhaustive: bool = True
    error: str = ""

    @property
    def warn(self) -> bool:
        limit = HAIRPIN_WARN_DG if self.kind == "hairpin" else DIMER_WARN_DG
        return self.found and self.dg <= limit


# --- parsing -----------------------------------------------------------------

_END_MARKS = ("5'-", "-3'", "5'", "3'")
# Invisible characters that str.isspace() does not cover (zero-width space/joiners, BOM).
_INVISIBLE = dict.fromkeys(map(ord, "​‌‍⁠﻿"), None)


def _strip_end_marks(text: str) -> str:
    text = text.strip()
    for mark in _END_MARKS:
        text = text.removeprefix(mark) if mark.startswith("5") else text.removesuffix(mark)
    return text.strip()


def _is_seq_token(token: str) -> bool:
    token = _strip_end_marks(token).upper()
    return bool(token) and all(c in IUPAC or c in "+U" for c in token)


def _char_name(c: str) -> str:
    return c if c.isascii() and c.isprintable() else f"U+{ord(c):04X}"


def _split_line(line: str) -> tuple[str | None, str]:
    """Split a stripped input line into (label, sequence); sequence "" = name-only line.

    Any whitespace inside the sequence (spaces, tabs, non-breaking spaces) is ignored later,
    so codon-style groups like 'GGT GGT TTC' or 'GGT<TAB>GGT' are one sequence.
    - Only sequence-like groups: a bare sequence.
    - ':', ';' or ',' or a tab separate label and sequence ('LNA probe<TAB>TG+A TTC').
    - Otherwise the first word is the label and the rest the sequence ('P1 GGT GGT TTC').
      If the rest is not sequence-like, the whole line is a label ('LNA probe').
    """
    tokens = line.split()
    if all(_is_seq_token(t) for t in tokens):
        return None, line
    cuts = [i for i in (line.find(sep) for sep in ":;,\t") if i >= 0]
    if cuts:  # the label ends at the first separator
        return line[:min(cuts)].strip(), line[min(cuts) + 1:]
    rest = tokens[1:]
    if rest and not all(_is_seq_token(t) for t in rest):
        return line, ""
    return tokens[0], " ".join(rest)


def parse_sequence(raw: str, label: str) -> tuple[str, frozenset[int]]:
    """Plain IUPAC sequence and LNA positions; whitespace is ignored, '+N' marks an LNA base."""
    seq: list[str] = []
    lna: set[int] = set()
    pending_lna = False
    bad: set[str] = set()
    for c in _strip_end_marks(raw.translate(_INVISIBLE)):
        if c.isspace():
            continue
        if c == "+":
            if pending_lna:
                raise ValueError(f"{label}: '++' is not valid LNA notation")
            pending_lna = True
            continue
        base = c.upper().replace("U", "T")
        if base not in IUPAC:
            bad.add(_char_name(c))
            continue
        c = base
        if pending_lna:
            lna.add(len(seq))
            pending_lna = False
        seq.append(c)
    if bad:
        raise ValueError(f"{label}: invalid characters {' '.join(sorted(bad))}")
    if pending_lna:
        raise ValueError(f"{label}: '+' at the end is not followed by a base")
    if not seq:
        raise ValueError(f"{label}: empty sequence")
    return "".join(seq), frozenset(lna)


def _line_entries(text: str) -> list[tuple[str | None, str]]:
    """(label, sequence) from line input; a name-only line labels the sequence on the next line."""
    entries: list[tuple[str | None, str]] = []
    pending: tuple[int, str] | None = None  # (line number, name) waiting for its sequence

    def missing(line_no: int, name: str) -> ValueError:
        return ValueError(ascii_only(
            f"line {line_no}: '{name}' has no sequence. Put the sequence after the name on "
            "the same line, or alone on the next line (if this line should contain a sequence, "
            "check it for invalid characters)."))

    for line_no, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()  # also removes trailing tabs and non-breaking spaces
        if not line or line.startswith("#"):
            continue
        name, seq = _split_line(line)
        if not seq.strip():  # name-only line
            if pending:
                raise missing(*pending)
            pending = (line_no, name or "")
        elif pending:
            if name is not None:
                raise missing(*pending)
            entries.append((pending[1], seq))
            pending = None
        else:
            entries.append((name, seq))
    if pending:
        raise missing(*pending)
    return entries


def _normalize_fasta(text: str) -> str:
    """Strip headers and remove all whitespace inside sequence lines before Bio.SeqIO.

    Bio.Seq only accepts ASCII, so a non-breaking space or other non-ASCII character in a
    sequence line is removed (whitespace) or reported here with the record name.
    """
    out, title = [], None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            title = line[1:].strip()
            out.append(line)
        elif title is not None and not line.startswith(";"):
            line = "".join(line.split())
            bad = sorted({_char_name(c) for c in line if not c.isascii()})
            if bad:
                raise ValueError(ascii_only(f"{title}: invalid characters {' '.join(bad)}"))
            out.append(line)
        else:
            out.append(line)
    return "\n".join(out)


def parse_sequences(text: str) -> list[Oligo]:
    """Parse FASTA (via Bio.SeqIO) or line input ('label<tab|:|;|,|space>sequence', bare
    sequences, or a name line followed by a sequence line)."""
    text = text.translate(_INVISIBLE)
    if any(line.lstrip().startswith(">") for line in text.splitlines()):
        # 'fasta-pearson' allows comment lines before the first record and ';' lines.
        # ('fasta-blast' in Biopython 1.88 repeats the first title for every record.)
        records = SeqIO.parse(io.StringIO(_normalize_fasta(text)), "fasta-pearson")
        entries = [(rec.description.strip() or None, str(rec.seq)) for rec in records]
    else:
        entries = _line_entries(text)

    oligos = []
    for number, (name, raw_seq) in enumerate(entries, 1):
        label = ascii_only(name or f"seq{number}")
        seq, lna = parse_sequence(raw_seq, f"#{number} {label}")
        oligos.append(Oligo(number, label, seq, lna))
    return oligos


# --- degenerate handling -----------------------------------------------------

def variant_combos(seqs: list[str], cap: int) -> tuple[Iterator[tuple[str, ...]], int, int, bool]:
    """Concrete variant combinations of degenerate sequences.

    Returns (iterator, total, checked, exhaustive): all combinations when there are at most
    `cap`, otherwise a reproducible random sample of `cap` combinations.
    """
    options = [[IUPAC[c] for c in s] for s in seqs]
    total = math.prod(len(o) for opts in options for o in opts)
    if total <= cap:
        per_seq = (map("".join, itertools.product(*opts)) for opts in options)
        return itertools.product(*per_seq), total, total, True
    rng = random.Random("|".join(seqs))
    sample = (tuple("".join(rng.choice(o) for o in opts) for opts in options) for _ in range(cap))
    return sample, total, cap, False


# --- Tm, GC, length, reverse complement --------------------------------------

def gc_content(seq: str) -> MinMeanMax:
    """Min, mean and max GC % over all variants (exact, no enumeration needed)."""
    gc = [sum(b in "GC" for b in IUPAC[c]) / len(IUPAC[c]) for c in seq]
    n = len(seq)
    return MinMeanMax(100 * sum(f == 1 for f in gc) / n, 100 * sum(gc) / n,
                      100 * sum(f > 0 for f in gc) / n)


def tm_biopython(seq: str, cond: Conditions) -> float:
    """Bio.SeqUtils.MeltingTemp.Tm_NN with the same salt and oligo concentrations.

    dnac1 = dnac2 = oligo/2 gives the same concentration term (C/4) as primer3.
    """
    return MeltingTemp.Tm_NN(seq, Na=cond.mv_conc, Mg=cond.dv_conc, dNTPs=cond.dntp_conc,
                             dnac1=cond.dna_conc / 2, dnac2=cond.dna_conc / 2,
                             nn_table=BIO_NN_TABLE, saltcorr=BIO_SALTCORR)


class _Stats:
    def __init__(self) -> None:
        self.lo, self.hi, self.total, self.n = math.inf, -math.inf, 0.0, 0

    def add(self, x: float) -> None:
        self.lo, self.hi = min(self.lo, x), max(self.hi, x)
        self.total, self.n = self.total + x, self.n + 1

    def result(self) -> MinMeanMax:
        return MinMeanMax(self.lo, self.total / self.n, self.hi)


def analyze(oligo: Oligo, cond: Conditions) -> BasicResult:
    """Tm (primer3, Biopython, LNA-corrected), GC %, reverse complement over all variants.

    The LNA-corrected Tm is primer3's DNA Tm plus the shift that the LNA nearest-neighbor
    parameters cause in the same model (identical to the full LNA model up to 60 nt).
    """
    combos, _, checked, exhaustive = variant_combos([oligo.seq], MAX_TM_VARIANTS)
    p3, bio, lna = _Stats(), _Stats(), _Stats()
    kw = cond.tm_kwargs()
    for (v,) in combos:
        tm = primer3.calc_tm(v, **kw)
        p3.add(tm)
        bio.add(tm_biopython(v, cond))
        if oligo.lna:
            lna.add(tm + tm_nn(v, oligo.lna, **kw) - tm_nn(v, frozenset(), **kw))
    return BasicResult(
        oligo=oligo,
        tm_primer3=p3.result(), tm_bio=bio.result(),
        tm_lna=lna.result() if oligo.lna else None,
        tm_checked=checked, tm_exhaustive=exhaustive,
        gc=gc_content(oligo.seq),
        revcomp=reverse_complement(oligo.seq),  # Bio.Seq handles IUPAC codes
    )


# --- hairpins and dimers (plain DNA; primer3 cannot model LNA) ---------------

def _most_stable(oligos: tuple[Oligo, ...], kind: str, seqs: list[str], calc, join: bool,
                 cap: int) -> StructResult:
    """Run `calc` on every variant combination and keep the most stable structure."""
    combos, total, checked, exhaustive = variant_combos(seqs, cap)
    res = StructResult(oligos, kind, n_total=total, n_checked=checked, exhaustive=exhaustive)
    best: tuple[float, tuple[str, ...]] | None = None
    try:
        for combo in combos:
            args = ("".join(combo),) if join else combo
            r = calc(*args, output_structure=False)
            if r.structure_found and (best is None or r.dg < best[0]):
                best = (r.dg, combo)
        if best:
            combo = best[1]
            r = calc(*(("".join(combo),) if join else combo), output_structure=True)
            res.found, res.tm, res.variant = True, r.tm, combo
            res.dg, res.dh, res.ds = r.dg / 1000, r.dh / 1000, r.ds
            res.structure = [ln.split("\t", 1)[-1] for ln in r.ascii_structure_lines
                             if ln.strip()]
    except (RuntimeError, ValueError) as exc:
        msg = exc.args[0] if exc.args else exc
        res.error = msg.decode() if isinstance(msg, bytes) else str(msg)
    return res


def hairpin(oligos: tuple[Oligo, ...], cond: Conditions) -> StructResult:
    """Hairpin of one oligo, or of several joined 5'->3'."""
    seqs = [o.seq for o in oligos]
    if sum(map(len, seqs)) > THAL_MAX_LEN:
        return StructResult(oligos, "hairpin",
                            error=f"longer than {THAL_MAX_LEN} nt (primer3 limit)")
    calc = lambda s, **kw: primer3.calc_hairpin(s, **cond.thermo_kwargs(), **kw)
    return _most_stable(oligos, "hairpin", seqs, calc, True, MAX_HAIRPIN_VARIANTS)


def hairpins(oligos: list[Oligo], cond: Conditions) -> tuple[list[StructResult], list[StructResult]]:
    """Hairpins of each oligo and of every ordered pair joined 5'->3' (A+B and B+A)."""
    single = [hairpin((o,), cond) for o in oligos]
    joined = [hairpin(pair, cond) for pair in itertools.permutations(oligos, 2)]
    return single, joined


def dimers(oligos: list[Oligo], cond: Conditions) -> tuple[list[StructResult], list[StructResult]]:
    """Homodimer of every oligo and heterodimer of every unordered pair."""
    homo_calc = lambda s, **kw: primer3.calc_homodimer(s, **cond.thermo_kwargs(), **kw)
    hetero_calc = lambda a, b, **kw: primer3.calc_heterodimer(a, b, **cond.thermo_kwargs(), **kw)
    homo = [_most_stable((o,), "dimer", [o.seq], homo_calc, False, MAX_DIMER_VARIANTS)
            for o in oligos]
    hetero = [_most_stable((a, b), "dimer", [a.seq, b.seq], hetero_calc, False,
                           MAX_DIMER_VARIANTS)
              for a, b in itertools.combinations(oligos, 2)]
    return homo, hetero
