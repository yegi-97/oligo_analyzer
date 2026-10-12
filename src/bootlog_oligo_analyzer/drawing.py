"""Draw primer3 hairpin and dimer structures on a Tkinter Canvas.

Each base is a labelled circle; consecutive bases are joined by a backbone line and paired
bases by a rung. Layout (pure geometry) is separated from rendering so it can be tested.
"""

import math
import tkinter as tk
from dataclasses import dataclass, field

D = 30  # distance between neighboring bases (px)
R = 11  # base circle radius (px)
RUNG_GAP = 2.2 * D  # distance between paired strands

FILL = {"A": "#a5d6a7", "T": "#ef9a9a", "G": "#ffe082", "C": "#90caf9"}  # by base (hairpins)
STRAND_COLORS = ["#1565c0", "#e65100"]  # first / second sequence: outline, backbone, labels
STRAND_FILLS = ["#bbdefb", "#ffe0b2"]  # by strand (dimers)
STRAND_NAMES = ["blue", "orange"]
LNA_RING = "#8e24aa"


@dataclass
class Node:
    x: float
    y: float
    base: str
    strand: int = 0  # 0 = first sequence, 1 = second sequence
    lna: bool = False


@dataclass
class Label:
    x: float
    y: float
    text: str
    strand: int
    anchor: str = "center"


@dataclass
class Layout:
    nodes: list[Node] = field(default_factory=list)
    backbone: list[tuple[int, int]] = field(default_factory=list)  # node index pairs
    rungs: list[tuple[int, int]] = field(default_factory=list)
    labels: list[Label] = field(default_factory=list)
    fill_by_strand: bool = False  # dimers: colour circles by strand instead of base


# --- hairpin -------------------------------------------------------------------

def hairpin_pairs(symbols: str) -> list[tuple[int, int]]:
    """Base pairs (i, j) from primer3's hairpin line: '/' = 5' side, '\\' = 3' side."""
    stack, pairs = [], []
    for k, c in enumerate(symbols):
        if c == "/":
            stack.append(k)
        elif c == "\\" and stack:
            pairs.append((stack.pop(), k))
    return sorted(pairs)


def _domains(pairs: list[tuple[int, int]]) -> list[list[tuple[int, int]]]:
    """Group pairs into stem-loops: each top-level pair with the pairs nested inside it."""
    domains: list[list[tuple[int, int]]] = []
    for i, j in sorted(pairs):
        if domains and i < domains[-1][0][1]:
            domains[-1].append((i, j))
        else:
            domains.append([(i, j)])
    return domains


def _stem_loop(dpairs: list[tuple[int, int]], x0: float) -> dict[int, tuple[float, float]] | None:
    """One stem-loop with its outermost pair at y = 0 and its 5' strand at x = x0.

    Stem rows go up from the outermost to the innermost pair; the loop bases sit on a
    semicircle on top. Returns None if the pairs are not one nested chain.
    """
    partner = {i: j for i, j in dpairs} | {j: i for i, j in dpairs}
    (i0, j0), (ii, jj) = dpairs[0], dpairs[-1]
    loop = list(range(ii + 1, jj))
    r = max(D, (len(loop) + 1) * D / math.pi)  # loop bases evenly on a semicircle
    width = 2 * r
    pos: dict[int, tuple[float, float]] = {}
    row, i, j = 0, i0, j0
    while i < j:
        if partner.get(i) == j:
            pos[i], pos[j] = (x0, -row * D), (x0 + width, -row * D)
            if (i, j) == (ii, jj):
                break
            i, j = i + 1, j - 1
        elif i not in partner and j not in partner:  # internal loop: unpaired, same row
            pos[i], pos[j] = (x0, -row * D), (x0 + width, -row * D)
            i, j = i + 1, j - 1
        elif i not in partner:  # bulge on the 5' strand
            pos[i] = (x0, -row * D)
            i += 1
        elif j not in partner:  # bulge on the 3' strand
            pos[j] = (x0 + width, -row * D)
            j -= 1
        else:  # both paired, but not with each other: branched, not a single stem
            return None
        row += 1
    else:
        return None
    top = -row * D
    for k, idx in enumerate(loop, 1):
        angle = math.pi - k * math.pi / (len(loop) + 1)
        pos[idx] = (x0 + r + r * math.cos(angle), top - r * math.sin(angle))
    return pos


def hairpin_layout(seq: str, symbols: str, lna: frozenset[int] = frozenset(),
                   boundary: int | None = None) -> Layout:
    """Stem as two vertical strands (5' left going up, 3' right coming down), loop on top.

    One stem-loop: unpaired ends hang down below the stem. Several stem-loops (primer3 can
    predict more than one): they stand side by side on a baseline, joined by their linker
    bases. `boundary` = length of the first sequence when two sequences are joined.
    """
    n = len(seq)
    pairs = hairpin_pairs(symbols)
    domains = _domains(pairs)
    pos: dict[int, tuple[float, float]] | None = {}

    if len(domains) == 1:
        pos = _stem_loop(domains[0], 0.0)
        if pos is not None:
            (i0, j0) = domains[0][0]
            width = pos[j0][0]
            for t, idx in enumerate(range(i0 - 1, -1, -1), 1):  # 5' tail hangs down left
                pos[idx] = (0.0, t * D)
            for t, idx in enumerate(range(j0 + 1, n), 1):  # 3' tail hangs down right
                pos[idx] = (width, t * D)
    elif domains:
        x = 0.0
        first = domains[0][0][0]
        for idx in range(first):  # 5' tail on the baseline, left of the first stem
            pos[idx] = (x - (first - idx) * D, 0.0)
        for k, dom in enumerate(domains):
            part = _stem_loop(dom, x)
            if part is None:
                pos = None
                break
            pos |= part
            end = dom[0][1]
            nxt = domains[k + 1][0][0] if k + 1 < len(domains) else n
            x_right = part[end][0]
            for t, idx in enumerate(range(end + 1, nxt), 1):  # linker / 3' tail on baseline
                pos[idx] = (x_right + t * D, 0.0)
            x = x_right + (nxt - end) * D

    if not pairs or pos is None or len(pos) != n:  # no stem or unexpected: plain strand
        pairs, pos = [], {i: (i * D, 0.0) for i in range(n)}

    lay = Layout()
    for i in range(n):
        strand = 1 if boundary is not None and i >= boundary else 0
        lay.nodes.append(Node(*pos[i], seq[i], strand, i in lna))
    lay.backbone = [(i, i + 1) for i in range(n - 1)]
    lay.rungs = pairs
    lay.labels = _end_labels(lay.nodes[0], lay.nodes[-1], below=len(domains) == 1 and bool(pairs))
    return lay


def _end_labels(first: Node, last: Node, below: bool) -> list[Label]:
    if below:
        return [Label(first.x, first.y + D, "5'", first.strand),
                Label(last.x, last.y + D, "3'", last.strand)]
    return [Label(first.x - D, first.y, "5'", first.strand),
            Label(last.x + D, last.y, "3'", last.strand)]


# --- dimer ---------------------------------------------------------------------

def dimer_layout(lines: list[str], top_name: str = "", bottom_name: str = "",
                 top_lna: frozenset[int] = frozenset(), bottom_lna: frozenset[int] = frozenset(),
                 bottom_len: int = 0) -> Layout:
    """Two antiparallel strands from primer3's 4 dimer lines.

    Lines: top unpaired, top paired (5'->3'), bottom paired, bottom unpaired (3'->5').
    Paired bases share a column and are joined by rungs; unpaired bases sit slightly outward.
    Each strand has its own colour, its name at its 5' start and its own 5'/3' marks.
    """
    width = max(map(len, lines))
    l1, l2, l3, l4 = (ln.ljust(width) for ln in (lines + [""] * 4)[:4])
    lay = Layout(fill_by_strand=True)
    top: list[int] = []
    bottom: list[int] = []
    for c in range(width):
        x = c * D
        if l2[c].isalpha() or l1[c].isalpha():
            paired = l2[c].isalpha()
            base = l2[c] if paired else l1[c]
            lay.nodes.append(Node(x, 0.0 if paired else -0.6 * D, base, 0, len(top) in top_lna))
            top.append(len(lay.nodes) - 1)
        if l3[c].isalpha() or l4[c].isalpha():
            paired = l3[c].isalpha()
            base = l3[c] if paired else l4[c]
            # bottom strand runs 3'->5' left to right: k-th base = index len-1-k
            idx = bottom_len - 1 - len(bottom)
            lay.nodes.append(Node(x, RUNG_GAP if paired else RUNG_GAP + 0.6 * D, base, 1,
                                  idx in bottom_lna))
            bottom.append(len(lay.nodes) - 1)
        if l2[c].isalpha() and l3[c].isalpha():
            lay.rungs.append((len(lay.nodes) - 2, len(lay.nodes) - 1))
    lay.backbone = list(zip(top, top[1:])) + list(zip(bottom, bottom[1:]))
    gap = R + 6
    if top:  # top strand starts (5') on the left
        t0, t1 = lay.nodes[top[0]], lay.nodes[top[-1]]
        lay.labels += [Label(t0.x - gap, t0.y, f"{top_name}  5'".strip(), 0, "e"),
                       Label(t1.x + gap, t1.y, "3'", 0, "w")]
    if bottom:  # bottom strand starts (5') on the right
        b0, b1 = lay.nodes[bottom[0]], lay.nodes[bottom[-1]]
        lay.labels += [Label(b0.x - gap, b0.y, "3'", 1, "e"),
                       Label(b1.x + gap, b1.y, f"5'  {bottom_name}".strip(), 1, "w")]
    return lay


# --- rendering -------------------------------------------------------------------

def render(canvas: tk.Canvas, lay: Layout, margin: int = 40) -> None:
    """Draw a layout and fit the scroll region."""
    canvas.delete("all")
    nodes = lay.nodes
    for a, b in lay.backbone:
        canvas.create_line(nodes[a].x, nodes[a].y, nodes[b].x, nodes[b].y,
                           fill=STRAND_COLORS[nodes[a].strand], width=2)
    for a, b in lay.rungs:
        na, nb = nodes[a], nodes[b]
        dist = math.hypot(nb.x - na.x, nb.y - na.y) or 1
        ux, uy = (nb.x - na.x) / dist, (nb.y - na.y) / dist
        canvas.create_line(na.x + ux * R, na.y + uy * R, nb.x - ux * R, nb.y - uy * R,
                           fill="#455a64", width=3)
    for node in nodes:
        if node.lna:
            canvas.create_oval(node.x - R - 4, node.y - R - 4, node.x + R + 4, node.y + R + 4,
                               outline=LNA_RING, width=2)
        fill = (STRAND_FILLS[node.strand] if lay.fill_by_strand
                else FILL.get(node.base, "#e0e0e0"))
        canvas.create_oval(node.x - R, node.y - R, node.x + R, node.y + R,
                           fill=fill, outline=STRAND_COLORS[node.strand], width=2)
        canvas.create_text(node.x, node.y, text=node.base, font=("Consolas", 10, "bold"))
    for lb in lay.labels:
        canvas.create_text(lb.x, lb.y, text=lb.text, anchor=lb.anchor,
                           font=("Segoe UI", 10, "bold"), fill=STRAND_COLORS[lb.strand])
    x0, y0, x1, y1 = canvas.bbox("all") or (0, 0, 0, 0)
    canvas.configure(scrollregion=(x0 - margin, y0 - margin, x1 + margin, y1 + margin))
    canvas.xview_moveto(0)
    canvas.yview_moveto(0)
