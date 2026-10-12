"""Tkinter front end: sequence input, conditions, one tab per analysis, detail windows."""

import queue
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText

from bootlog_oligo_analyzer import core, drawing, info

EXAMPLE = """fwd\tACGTGGATCCAAGCTTGCAGT
rev\tTGCAAGCTTGGATCCACGTCA
degenerate\tGAYCARGTNTGGAARGC
LNA probe\tTG+A TTC +TCA G+CC C+TT C+GC
"""

INPUT_HELP = ("'name<TAB>sequence', 'name: sequence', 'name sequence', a name line followed by "
              "a sequence line, or FASTA. Spaces/tabs inside a sequence are ignored "
              "(GGT GGT TTC), +N = LNA, IUPAC codes allowed.")

CONDITIONS = [  # (attribute, label)
    ("mv_conc", "Na+/K+ mM"),
    ("dv_conc", "Mg2+ mM"),
    ("dntp_conc", "dNTP mM"),
    ("dna_conc", "Oligo nM"),
    ("temp_c", "Temp C"),
]

PROPERTY_COLUMNS = [  # (id, heading, width, anchor)
    ("num", "#", 45, "e"), ("name", "Name", 170, "w"), ("len", "Len", 50, "e"),
    ("lna", "LNA", 45, "e"), ("variants", "Variants", 80, "e"),
    ("tm_p3", "Tm primer3 (C)", 140, "e"), ("tm_bio", "Tm Biopython (C)", 140, "e"),
    ("tm_lna", "Tm LNA (C)", 140, "e"), ("gc", "GC (%)", 130, "e"),
]
STRUCT_COLUMNS = [
    ("num", "#", 80, "w"), ("name", "Name", 260, "w"), ("dg", "dG (kcal/mol)", 110, "e"),
    ("tm", "Tm (C)", 80, "e"), ("flag", "Flag", 120, "w"),
]


# --- formatting helpers ---------------------------------------------------------

def _range_cell(v: core.MinMeanMax | None, degenerate: bool) -> str:
    if v is None:
        return "-"
    return f"{v.mean:.1f} ({v.min:.1f}-{v.max:.1f})" if degenerate else f"{v.mean:.1f}"


def _ids(res: core.StructResult) -> tuple[str, str]:
    sep = " + " if res.kind == "hairpin" else " x "
    return (sep.join(f"#{o.number}" for o in res.oligos), sep.join(o.name for o in res.oligos))


def _flag(res: core.StructResult) -> tuple[str, str]:
    """(cell text, row tag)"""
    if res.error:
        return "not calculated", "dim"
    if not res.found:
        return "no structure", "ok"
    return ("WARN", "warn") if res.warn else ("ok", "")


def _variants_text(total: int, checked: int, exhaustive: bool) -> str:
    if total == 1:
        return "1 (not degenerate)"
    if exhaustive:
        return f"{total:,} (all evaluated)"
    return f"{total:,} ({checked:,} randomly sampled)"


def _strand_oligos(res: core.StructResult) -> list[core.Oligo]:
    """Oligo of each drawn strand (a homodimer has two strands of the same oligo)."""
    if res.kind == "dimer":
        return [res.oligos[0], res.oligos[-1]]
    return list(res.oligos)


def _struct_layout(res: core.StructResult) -> drawing.Layout:
    if res.kind == "hairpin":
        lna, offset = set(), 0
        for o in res.oligos:
            lna |= {offset + i for i in o.lna}
            offset += len(o.seq)
        boundary = len(res.oligos[0].seq) if len(res.oligos) > 1 else None
        symbols, seq = res.structure[0], res.structure[1]
        return drawing.hairpin_layout(seq, symbols, frozenset(lna), boundary)
    top, bottom = _strand_oligos(res)
    return drawing.dimer_layout(res.structure, top.label, bottom.label,
                                top.lna, bottom.lna, len(bottom.seq))


Lines = list[tuple[str, int | None]]  # (text, strand index for colouring or None)


def _struct_detail_text(res: core.StructResult) -> Lines:
    num, names = _ids(res)
    title = ("Hairpin" if res.kind == "hairpin" else
             "Homodimer" if len(res.oligos) == 1 else "Heterodimer")
    strands = _strand_oligos(res)
    width = max(len(o.label) for o in strands)
    lines: Lines = [(f"{title}: {num} {names}", None)]
    for k, o in enumerate(strands):
        lines.append((f"  {o.label:<{width}}  strand {k + 1} ({drawing.STRAND_NAMES[k]}):"
                      f"  5'-{o.marked}-3'", k))
    if res.error:
        lines.append((f"Not calculated: {res.error}", None))
    elif not res.found:
        lines.append(("No structure predicted.", None))
    else:
        limit = core.HAIRPIN_WARN_DG if res.kind == "hairpin" else core.DIMER_WARN_DG
        lines += [
            (f"dG {res.dg:.2f} kcal/mol   dH {res.dh:.2f} kcal/mol   "
             f"dS {res.ds:.1f} cal/(mol K)   Tm {res.tm:.1f} C", None),
            (f"Flag: {'WARN' if res.warn else 'ok'} (WARN if dG <= {limit} kcal/mol)", None),
        ]
    lines.append(("Variants checked: "
                  + _variants_text(res.n_total, res.n_checked, res.exhaustive), None))
    if res.found and res.n_total > 1:
        lines.append(("Most stable variant: " + " / ".join(res.variant), None))
    if any(o.lna for o in res.oligos):
        lines.append(("LNA bases (purple ring) are treated as DNA: primer3 cannot model LNA.",
                      None))
    if res.found and res.kind == "dimer":
        # primer3 alignment: lines 1-2 = first strand 5'->3', lines 3-4 = second strand 3'<-5'
        lines.append(("Alignment (paired bases in the middle two lines):", None))
        for k, ln in enumerate(res.structure[:4]):
            strand = 0 if k < 2 else 1
            o = strands[strand]
            ends = ("5'", "3'") if strand == 0 else ("3'", "5'")
            lines.append((f"  {o.label:<{width}}  {ends[0]}  {ln.rstrip()}", strand))
    return lines


def _property_detail_text(r: core.BasicResult, hp: core.StructResult) -> str:
    o = r.oligo
    rows = [("Tm primer3 (C)", r.tm_primer3), ("Tm Biopython (C)", r.tm_bio),
            ("Tm LNA (C)", r.tm_lna), ("GC (%)", r.gc)]
    lines = [
        o.label,
        f"Sequence 5'->3'      {o.marked}",
        f"Length               {len(o.seq)} nt",
        "LNA positions        " + (", ".join(str(i + 1) for i in sorted(o.lna)) or "-"),
        "Variants             " + _variants_text(o.n_variants, r.tm_checked, r.tm_exhaustive),
        "",
        f"{'':<20} {'min':>7} {'mean':>7} {'max':>7}",
    ]
    for label, v in rows:
        lines.append(f"{label:<20} " + (f"{v.min:>7.1f} {v.mean:>7.1f} {v.max:>7.1f}" if v
                                        else f"{'-':>7} {'-':>7} {'-':>7}"))
    lines += ["", f"Reverse complement   {r.revcomp}", ""]
    if hp.error:
        lines.append(f"Hairpin: not calculated ({hp.error})")
    elif not hp.found:
        lines.append("Hairpin: none predicted")
    else:
        lines.append(f"Hairpin: dG {hp.dg:.2f} kcal/mol, Tm {hp.tm:.1f} C"
                     + (f", most stable variant {hp.variant[0]}" if hp.n_total > 1 else "")
                     + (" (WARN)" if hp.warn else ""))
    return "\n".join(lines)


# --- windows ----------------------------------------------------------------------

def text_window(parent: tk.Misc, title: str, text: str, size: str = "640x480") -> None:
    win = tk.Toplevel(parent)
    win.title(title)
    win.geometry(size)
    box = ScrolledText(win, wrap="word", font=("Segoe UI", 10), padx=10, pady=8)
    box.pack(fill="both", expand=True)
    box.insert("1.0", text)
    box.configure(state="disabled")


def detail_window(parent: tk.Misc, title: str, lines: Lines,
                  layout: drawing.Layout | None) -> None:
    """Values as selectable text on top (strand lines in strand colour), drawing below."""
    win = tk.Toplevel(parent)
    win.title(title)
    win.geometry("900x700")
    box = tk.Text(win, height=min(len(lines), 20), wrap="none", font=("Consolas", 10),
                  padx=8, pady=6)
    box.pack(fill="x")
    for k, color in enumerate(drawing.STRAND_COLORS):
        box.tag_configure(f"strand{k}", foreground=color)
    for text, strand in lines:
        box.insert("end", text + "\n", (f"strand{strand}",) if strand is not None else ())
    box.configure(state="disabled")
    if layout is None:
        return
    frame = ttk.Frame(win)
    frame.pack(fill="both", expand=True)
    frame.rowconfigure(0, weight=1)
    frame.columnconfigure(0, weight=1)
    canvas = tk.Canvas(frame, background="white", highlightthickness=0)
    canvas.grid(row=0, column=0, sticky="nsew")
    vsb = ttk.Scrollbar(frame, orient="vertical", command=canvas.yview)
    hsb = ttk.Scrollbar(frame, orient="horizontal", command=canvas.xview)
    vsb.grid(row=0, column=1, sticky="ns")
    hsb.grid(row=1, column=0, sticky="ew")
    canvas.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
    drawing.render(canvas, layout)


class ResultTab(ttk.Frame):
    """A tab with a summary line, an Info button and a clickable result table."""

    def __init__(self, master, columns, info_text, on_open) -> None:
        super().__init__(master, padding=4)
        self.info_text, self.on_open = info_text, on_open
        self.cond: core.Conditions | None = None
        self.payloads: dict[str, object] = {}

        bar = ttk.Frame(self)
        bar.pack(fill="x", pady=(0, 4))
        self.summary = ttk.Label(bar, text="Not calculated yet.")
        self.summary.pack(side="left")
        ttk.Button(bar, text="Info", command=self.show_info).pack(side="right")

        frame = ttk.Frame(self)
        frame.pack(fill="both", expand=True)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        self.tree = ttk.Treeview(frame, columns=[c[0] for c in columns], show="headings",
                                 selectmode="browse", cursor="hand2")
        for cid, heading, width, anchor in columns:
            self.tree.heading(cid, text=heading, anchor=anchor)
            self.tree.column(cid, width=width, anchor=anchor, stretch=cid == "name")
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb = ttk.Scrollbar(frame, orient="vertical", command=self.tree.yview)
        vsb.grid(row=0, column=1, sticky="ns")
        self.tree.configure(yscrollcommand=vsb.set)
        self.tree.tag_configure("warn", foreground="#c62828")
        self.tree.tag_configure("ok", foreground="#2e7d32")
        self.tree.tag_configure("dim", foreground="#777777")
        self.tree.tag_configure("group", font=("Segoe UI", 9, "bold"), background="#eceff1")
        self.tree.bind("<ButtonRelease-1>", self._clicked)
        self.tree.bind("<Return>", lambda e: self._open(self.tree.focus()))

    def clear(self, summary: str) -> None:
        self.tree.delete(*self.tree.get_children())
        self.payloads.clear()
        self.summary.configure(text=summary)

    def add_group(self, title: str) -> None:
        self.tree.insert("", "end", values=("", title), tags=("group",))

    def add_row(self, values, payload, tag: str = "") -> None:
        iid = self.tree.insert("", "end", values=values, tags=(tag,) if tag else ())
        self.payloads[iid] = payload

    def _clicked(self, event) -> None:
        if self.tree.identify_region(event.x, event.y) == "cell":
            self._open(self.tree.identify_row(event.y))

    def _open(self, iid: str) -> None:
        if iid in self.payloads:
            self.on_open(self.payloads[iid], self.cond)

    def show_info(self) -> None:
        cond = self.cond or self.winfo_toplevel().current_conditions() or core.Conditions()
        text_window(self, "Info", self.info_text(cond))


class App(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("Oligo Analyzer")
        self.geometry("1200x800")
        self.minsize(800, 550)
        self._results: queue.Queue = queue.Queue()

        root = ttk.Frame(self, padding=8)
        root.pack(fill="both", expand=True)
        root.columnconfigure(0, weight=1)
        root.rowconfigure(4, weight=1)

        # Sequence input
        head = ttk.Frame(root)
        head.grid(row=0, column=0, sticky="ew")
        ttk.Label(head, text=INPUT_HELP).pack(side="left")
        ttk.Button(head, text="Clear", command=lambda: self.input.delete("1.0", "end")
                   ).pack(side="right")
        ttk.Button(head, text="Load file...", command=self.load_file).pack(side="right", padx=4)
        self.input = ScrolledText(root, height=7, font=("Consolas", 10), undo=True)
        self.input.grid(row=1, column=0, sticky="nsew", pady=(2, 6))
        self.input.insert("1.0", EXAMPLE)

        # Conditions
        cond_frame = ttk.LabelFrame(root, text="Conditions", padding=4)
        cond_frame.grid(row=2, column=0, sticky="ew")
        defaults = core.Conditions()
        self.cond_vars: dict[str, tk.StringVar] = {}
        for i, (attr, label) in enumerate(CONDITIONS):
            ttk.Label(cond_frame, text=label).grid(row=0, column=2 * i, padx=(8, 2))
            var = tk.StringVar(value=f"{getattr(defaults, attr):g}")
            ttk.Entry(cond_frame, textvariable=var, width=7).grid(row=0, column=2 * i + 1)
            self.cond_vars[attr] = var

        # Action buttons
        actions = ttk.Frame(root)
        actions.grid(row=3, column=0, sticky="ew", pady=6)
        self.buttons = []
        for text, tasks in [("Run all", ("props", "hairpins", "dimers")),
                            ("Properties", ("props",)),
                            ("Hairpins", ("hairpins",)),
                            ("Dimers", ("dimers",))]:
            btn = ttk.Button(actions, text=text, command=lambda t=tasks: self.start(t))
            btn.pack(side="left", padx=(0, 6), ipady=3)
            self.buttons.append(btn)
        self.status = ttk.Label(actions, foreground="#555555")
        self.status.pack(side="left", padx=8)

        # Result tabs
        self.notebook = ttk.Notebook(root)
        self.notebook.grid(row=4, column=0, sticky="nsew")
        self.tabs = {
            "props": ResultTab(self.notebook, PROPERTY_COLUMNS, info.properties,
                               self.open_property),
            "hairpins": ResultTab(self.notebook, STRUCT_COLUMNS, info.hairpins,
                                  self.open_structure),
            "dimers": ResultTab(self.notebook, STRUCT_COLUMNS, info.dimers,
                                self.open_structure),
        }
        for key, title in [("props", "Properties"), ("hairpins", "Hairpins"),
                           ("dimers", "Dimers")]:
            self.notebook.add(self.tabs[key], text=title)

    # --- input -------------------------------------------------------------------

    def load_file(self) -> None:
        path = filedialog.askopenfilename(
            filetypes=[("Sequence files", "*.fa *.fasta *.fas *.txt *.tsv *.csv"), ("All", "*")])
        if path:
            with open(path, encoding="utf-8", errors="replace") as fh:
                self.input.delete("1.0", "end")
                self.input.insert("1.0", fh.read())

    def current_conditions(self) -> core.Conditions | None:
        try:
            return core.Conditions(**{k: float(v.get().replace(",", "."))
                                      for k, v in self.cond_vars.items()})
        except ValueError:
            return None

    def read_input(self) -> tuple[list[core.Oligo], core.Conditions] | None:
        cond = self.current_conditions()
        if cond is None:
            messagebox.showerror("Invalid input", "Conditions must be numbers.")
            return None
        try:
            oligos = core.parse_sequences(self.input.get("1.0", "end"))
        except ValueError as exc:
            messagebox.showerror("Invalid input", str(exc))
            return None
        if not oligos:
            messagebox.showinfo("No sequences", "Please enter at least one sequence.")
            return None
        return oligos, cond

    # --- calculation (worker thread) ------------------------------------------------

    def start(self, tasks: tuple[str, ...]) -> None:
        parsed = self.read_input()
        if parsed is None:
            return
        oligos, cond = parsed
        self.notebook.select(self.tabs[tasks[0]])
        for t in tasks:
            self.tabs[t].clear("Calculating...")
        for b in self.buttons:
            b.state(["disabled"])
        self.status.configure(text="Calculating...")
        self.config(cursor="watch")

        def work() -> None:
            try:
                out = {}
                for t in tasks:
                    if t == "props":
                        out[t] = [core.analyze(o, cond) for o in oligos]
                    elif t == "hairpins":
                        out[t] = core.hairpins(oligos, cond)
                    else:
                        out[t] = core.dimers(oligos, cond)
                self._results.put((out, oligos, cond))
            except Exception as exc:  # report unexpected failures instead of hanging
                self._results.put(exc)

        threading.Thread(target=work, daemon=True).start()
        self.after(100, self._poll)

    def _poll(self) -> None:
        try:
            result = self._results.get_nowait()
        except queue.Empty:
            self.after(100, self._poll)
            return
        for b in self.buttons:
            b.state(["!disabled"])
        self.config(cursor="")
        self.status.configure(text="")
        if isinstance(result, Exception):
            messagebox.showerror("Calculation failed", f"{type(result).__name__}: {result}")
            return
        out, oligos, cond = result
        if "props" in out:
            self.fill_properties(out["props"], cond)
        if "hairpins" in out:
            single, joined = out["hairpins"]
            self.fill_structures("hairpins", cond, [("Individual sequences", single),
                                                    ("Joined pairs (A + B)", joined)])
        if "dimers" in out:
            homo, hetero = out["dimers"]
            self.fill_structures("dimers", cond, [("Homodimers", homo),
                                                  ("Heterodimers", hetero)])

    # --- tables ------------------------------------------------------------------------

    def fill_properties(self, results: list[core.BasicResult], cond: core.Conditions) -> None:
        tab = self.tabs["props"]
        tab.cond = cond
        n_lna = sum(bool(r.oligo.lna) for r in results)
        n_deg = sum(r.oligo.n_variants > 1 for r in results)
        tab.clear(f"{len(results)} sequences, {n_deg} degenerate, {n_lna} with LNA")
        for r in results:
            o = r.oligo
            deg = o.n_variants > 1
            tab.add_row((f"#{o.number}", o.name, len(o.seq), len(o.lna) or "-",
                         f"{o.n_variants:,}", _range_cell(r.tm_primer3, deg),
                         _range_cell(r.tm_bio, deg), _range_cell(r.tm_lna, deg),
                         _range_cell(r.gc, deg)), r)

    def fill_structures(self, key: str, cond: core.Conditions,
                        groups: list[tuple[str, list[core.StructResult]]]) -> None:
        tab = self.tabs[key]
        tab.cond = cond
        all_results = [r for _, results in groups for r in results]
        n_warn = sum(r.warn for r in all_results)
        tab.clear(f"{len(all_results)} predictions, {n_warn} WARN")
        for title, results in groups:
            tab.add_group(title)
            if not results:
                tab.tree.insert("", "end", values=("", "(needs at least two sequences)"),
                                tags=("dim",))
            for r in results:
                num, names = _ids(r)
                flag, tag = _flag(r)
                tab.add_row((num, names, f"{r.dg:.2f}" if r.found else "-",
                             f"{r.tm:.1f}" if r.found else "-", flag), r, tag)

    # --- detail views ---------------------------------------------------------------------

    def open_property(self, r: core.BasicResult, cond: core.Conditions) -> None:
        self.config(cursor="watch")
        self.update_idletasks()
        try:
            hp = core.hairpin((r.oligo,), cond)
        finally:
            self.config(cursor="")
        layout = _struct_layout(hp) if hp.found else None
        text = _property_detail_text(r, hp)
        detail_window(self, f"{r.oligo.label} - properties",
                      [(line, None) for line in text.splitlines()], layout)

    def open_structure(self, res: core.StructResult, cond: core.Conditions) -> None:
        num, names = _ids(res)
        layout = _struct_layout(res) if res.found else None
        detail_window(self, f"{num} {names} - {res.kind}", _struct_detail_text(res), layout)


def run() -> None:
    App().mainloop()
