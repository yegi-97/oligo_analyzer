# bootlog-oligo-analyzer

Slim Tkinter tool for oligo/primer analysis based on [primer3-py](https://github.com/libnano/primer3-py).

## Run

```
uv run bootlog-oligo-analyzer
```

## Input

Paste sequences straight from an ordering sheet or spreadsheet:

| Format | Example |
| --- | --- |
| tab-separated (sheet) | `LNA probe<TAB>TG+A TTC +TCA G+CC C+TT C+GC` |
| `name: sequence` (also `;` or `,`) | `probe 1: 5'-TG+A TTC +TCA-3'` |
| `name sequence` | `P1 GGT GGT TTC TGG GGT` |
| name line, then sequence line | `fwd<TAB><TAB>` followed by `GGT GGT TTC TGG GGT` |
| bare sequence | `GGT GGT TTC TGG GGT` (named `seqN`) |
| FASTA (Bio.SeqIO) | `>probe A` followed by sequence lines |

- All whitespace inside a sequence line is ignored: spaces (also several), tabs, non-breaking
  and zero-width spaces, so codon-style triplets work. Trailing tabs or spaces after a name
  (spreadsheet paste) are fine. A sequence must be on one line (except in FASTA).
- The name ends at the first tab, `:`, `;` or `,`, so names may contain spaces
  (`LNA probe<TAB>...`). Everything after it is sequence; extra spreadsheet columns after the
  sequence are reported as invalid characters, so copy only the name and sequence columns.
- `5'-`/`-3'` marks are stripped, `U` is read as `T`.
- `+N` marks an LNA base. IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) are allowed.
- Names are labels only. Every analysis runs on every sequence and every pair, in input order;
  nothing is ever matched, paired or grouped by name. Sequences are identified as `#n name`.
  Without a separator, a first word made only of base letters (e.g. `CAT`) is read as
  sequence, so use a tab or `:` for such names.

## Tabs

Each button fills and opens its own tab; **Run all** fills all three. Every tab has an
**Info** button with the column definitions, methods, references and the conditions used.
Click a row to open a detail window with all values and the structure drawing.

- **Properties**: length, LNA count, variants, Tm primer3, Tm Biopython, Tm LNA, GC %
  (degenerate sequences: mean with min-max range). Detail: min/mean/max table, LNA positions,
  reverse complement, predicted hairpin.
- **Hairpins**: each sequence and every joined ordered pair (A + B). Detail: dG, dH, dS, Tm,
  most stable variant, and the hairpin drawn as stem (two parallel strands with rungs) and a
  semicircular loop, each base a labelled circle, 5'/3' ends marked, LNA bases ringed. When
  primer3 predicts several stem-loops, they stand side by side on a baseline.
- **Dimers**: homodimers and heterodimers. Detail: values and the two strands drawn
  antiparallel with the paired region aligned. Each strand has its own colour (first sequence
  blue, second orange), its name at its 5' start and its own 5'/3' marks; in the detail text
  every strand line, including the alignment, starts with the sequence name in that colour.

## Models

- **Tm primer3**: `calc_tm` (SantaLucia 1998 nearest-neighbor, SantaLucia salt correction,
  Mg2+ converted to Na+ equivalents). LNA bases are counted as DNA.
- **Tm Biopython**: `Bio.SeqUtils.MeltingTemp.Tm_NN` with the same salt and oligo
  concentrations (`dnac1 = dnac2 = oligo/2`, the same C/4 term as primer3), table `DNA_NN3`
  and salt correction method 7 (Owczarzy 2008). With Biopython's default salt correction it is
  identical to primer3, so the column shows the effect of the Mg2+-aware salt model.
  Set `BIO_SALTCORR` / `BIO_NN_TABLE` in `core.py` to compare other models.
- **Tm LNA**: the primer3 model with every step that touches an LNA replaced by the LNA/DNA
  parameters of Owczarzy et al. 2011, *Biochemistry* 50:9352, Table 2 (`lna.py`). For plain DNA
  it reproduces primer3 exactly; above 60 nt the LNA shift is added to primer3's Tm.
- **Reverse complement**: `Bio.Seq.reverse_complement` (IUPAC-aware). FASTA: `Bio.SeqIO`
  (`fasta-pearson`; Biopython 1.88's `fasta-blast` repeats the first title for every record).
- **Hairpins/dimers**: primer3 on the plain DNA sequence, since primer3 cannot model LNA.
  Hairpins need <= 60 nt; for dimers at least one of the two sequences must be <= 60 nt.
- **Degenerate sequences**: all variants (Tm up to 50,000, hairpins up to 100,000, dimers up
  to 20,000 combinations; above that a reproducible random sample). Hairpins and dimers report
  the most stable variant. GC % min/mean/max is exact.

All generated text is ASCII (non-ASCII characters in names become `?`).

## Build a Windows executable

```
uv run pyinstaller OligoAnalyzer.spec
```

Produces `dist/OligoAnalyzer.exe` (single file, no console). The spec was generated with
`--onefile --windowed --collect-all primer3`; `--collect-all primer3` is required because
primer3 silently drops its functions when its compiled submodules are missing. Biopython pulls
in numpy, so the exe is about 27 MB.

## license

not implemented