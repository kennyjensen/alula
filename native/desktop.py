#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Alula desktop viewer for the native panel and Euler solvers."""
import json
import math
from pathlib import Path
import subprocess
import tempfile
import threading
try:
    import tkinter as tk
    from tkinter import filedialog, ttk
except ImportError:
    raise SystemExit("Install the Ubuntu package python3-tk to run the Alula desktop GUI.") from None

HERE = Path(__file__).resolve().parent
AIRFOILS = HERE.parent / "third_party" / "airfoils"
PRESETS = ("NACA 0012", "RAE 2822", "NLR 7301 + flap")
PANEL_MODE = "Incompressible panel"
VISCOUS_MODE = "Incompressible panel + BL"
PRESSURE_VIEW, BL_VIEW = "Pressure distribution", "Displacement thickness"
BG, PAPER, INK, MUTED = "#edf2f4", "#ffffff", "#192f3a", "#526974"
UPPER, LOWER, GRID, LIVE = "#05766f", "#b24e23", "#dde5e8", "#996018"
FIELDS = ("Alpha (°)", "Mach", "Reynolds number", "Ncrit", "Upper trip (x/c)",
          "Lower trip (x/c)", "Surface intervals", "Tubes per passage")
DEFAULTS = ((2, .2, 1000000, 9, 1, 1, 8, 4),
            (2.68, .74, 2700000, 4, .03, .07, 16, 4),
            (6, .185, 2510000, 9, 1, 1, 10, 4))


def preset_geometry(name):
    if name == PRESETS[1]:
        lines = (AIRFOILS / "rae2822-mses/blade.rae").read_text().splitlines()
        return [("S", [tuple(map(float, row.split())) for row in lines[2:] if row.strip()])]
    if name == PRESETS[2]:
        elements = json.loads((AIRFOILS / "nlr7301/geometry.json").read_text())["elements"]
        return [(f"F {e['points'].index(e['upperTrailingEdge'])} {e['points'].index(e['lowerTrailingEdge'])}",
                 e["points"]) for e in elements]
    upper = []
    for i in range(81):
        x = .5 * (1 + math.cos(math.pi * i / 80))
        y = .6 * (.2969 * math.sqrt(x) - .1260*x - .3516*x*x + .2843*x**3 - .1036*x**4)
        upper.append((x, y if 0 < i < 80 else 0.0))
    return [("S", upper + [(x, -y) for x, y in reversed(upper[1:-1])] + [upper[0]])]


def read_contours(path):
    words = iter(path.read_text().split())
    elements = []
    try:
        count = int(next(words))
        if not 1 <= count <= 6:
            raise ValueError
        for _ in range(count):
            size, token = int(next(words)), next(words)
            if size < 4:
                raise ValueError
            kind = token if token in ("S", "F") else "S"
            if token == "F":
                a, b = int(next(words)), int(next(words))
                if not (0 <= a < size and 0 <= b < size and a != b):
                    raise ValueError
                kind = f"F {a} {b}"
            first = next(words) if token in ("S", "F") else token
            points = [(float(first), float(next(words)))]
            for _ in range(size-1):
                points.append((float(next(words)), float(next(words))))
            if not all(math.isfinite(v) for p in points for v in p):
                raise ValueError
            elements.append((kind, points))
        if next(words, None) is not None:
            raise ValueError
    except (ValueError, StopIteration):
        raise ValueError("Use a native contour file with an element count, point counts and coordinate pairs.") from None
    return elements


def solver_args(name, mode, values, elements):
    if mode == PANEL_MODE:
        try:
            alpha = float(values[0])
        except ValueError:
            raise ValueError("Enter a numeric angle of attack.") from None
        if not -90 <= alpha <= 90:
            raise ValueError("Alpha must be within ±90°.")
        return [str(HERE / "alula-panel"), "--topology", "--cp", str(alpha)]
    if mode == VISCOUS_MODE:
        try:
            alpha, re, ncrit, top, bottom = map(float, (values[0], *values[2:6]))
        except ValueError:
            raise ValueError("Enter numeric flow conditions and trip fractions.") from None
        if not (-20 <= alpha <= 20 and 5e4 <= re <= 5e7 and 0 <= ncrit <= 14
                and 0 < top <= 1 and 0 < bottom <= 1):
            raise ValueError("Use α within ±20°, Re from 50,000 to 50 million, Ncrit 0–14 and trips in (0, 1].")
        return ([str(HERE / "alula-viscous"), "--topology", "--live-progress"]
                + ["--trip", str(top), str(bottom)] * len(elements) + list(map(str, (alpha, re, ncrit))))
    try:
        coupled = mode == "Euler + boundary layer"
        conditions = list(map(float, values[:6] if coupled else values[:2]))
        alpha, mach = conditions[:2]
        intervals, tubes = map(int, values[6:])
        if not all(math.isfinite(v) for v in conditions):
            raise ValueError
    except ValueError:
        raise ValueError("Enter finite numbers and whole-number mesh counts.") from None
    if not (-90 <= alpha <= 90 and 0 < mach < 1):
        raise ValueError("Alpha must be within ±90° and Mach must be between 0 and 1.")
    if not (4 <= intervals <= 512 and 1 <= tubes <= 256):
        raise ValueError("Use 4–512 surface intervals and 1–256 tubes.")
    args = [str(HERE / "alula-euler"), "--topology", "--cp", "--live-progress", "--potential-grid"]
    if intervals > 8:
        args += ["--coarse-start", "8x4"]
    fine_rae = name == PRESETS[1] and intervals >= 64
    if fine_rae:
        args += ["--domain", "-2", "3", "-3", "3.5", "--surface-uniformity", ".7"]
        if intervals == 64:
            args += ["--inlet-intervals", "16", "--wake-intervals", "16"]
    if coupled:
        re, ncrit, top, bottom = conditions[2:]
        if re <= 0 or ncrit < 0 or not (0 < top <= 1 and 0 < bottom <= 1):
            raise ValueError("Use positive Re, nonnegative Ncrit and trip fractions in (0, 1].")
        args += ["--coupled", str(re), str(ncrit)]
        args += ["--trip", str(top), str(bottom)] * len(elements)
        if fine_rae and mach > .64:
            args += ["--mach-start", ".64"]
        if name == PRESETS[2] and intervals > 8:
            args += ["--refine-coupled"]
    return args + list(map(str, (alpha, mach, intervals, tubes)))


def panel_pressure(elements, values):
    """Map the CLI's input-order panel averages to material surface curves."""
    if len(values) != sum(len(points)-1 for _, points in elements):
        raise ValueError("Panel pressure count differs from the input geometry.")
    samples, offset = [], 0
    for contour, (kind, points) in enumerate(elements):
        n = len(points)-1
        upper, lower = (0, n) if kind == "S" else map(int, kind.split()[1:])
        wetted = n if kind == "S" else (lower-upper) % n
        te = tuple((a+b)/2 for a, b in zip(points[upper], points[lower]))
        nose = max(range(wetted+1), key=lambda k: math.dist(points[(upper+k) % n], te))
        for k in range(n):
            i = (upper+k) % n
            samples.append({"contour": contour, "surface": "base" if k >= wetted else
                            ("upper" if k < nose else "lower"),
                            "x": (points[i][0] + points[i+1][0])/2, "cp": values[offset+i]})
        offset += n
    return samples


class Updates:
    """Keep every residual and just the most recent plot between GUI repaints."""
    def __init__(self):
        self.lock = threading.Lock()
        self.samples, self.latest, self.finished = [], None, None
        self.plot_data = {}

    def iterate(self, frame):
        with self.lock:
            self.samples.append((frame["stage"], frame["mach"], frame["iteration"], frame["residual"]))
            if frame.get("plotsDeferred"):
                frame.update(self.plot_data)
            else:
                self.plot_data = {k: frame[k] for k in ("cp", "streamlines", "boundaryLayerCurves") if k in frame}
            self.latest = frame

    def finish(self, code, report, detail):
        with self.lock:
            self.finished = (code, report, detail)

    def take(self):
        with self.lock:
            data = self.samples, self.latest, self.finished
            self.samples, self.latest, self.finished = [], None, None
            return data


def launch(args, elements, events):
    # The CLI streams accepted iterates on stderr and writes one final report
    # on stdout. A file for the report avoids any pipe-capacity dependency.
    with tempfile.TemporaryFile(mode="w+") as geometry:
        print(len(elements), file=geometry)
        for kind, points in elements:
            print(len(points), kind, file=geometry)
            for x, y in points:
                print(x, y, file=geometry)
        geometry.seek(0)
        report = tempfile.TemporaryFile(mode="w+")
        try:
            process = subprocess.Popen(args, stdin=geometry, stdout=report, stderr=subprocess.PIPE,
                                       text=True, encoding="utf-8", errors="replace")
        except OSError:
            report.close()
            raise

    def collect():
        detail = ""
        retained = {}
        with process.stderr, report:
            for line in process.stderr:
                try:
                    event = json.loads(line)
                    if event.get("type") == "iterate":
                        events.iterate(event)
                    elif event.get("type") == "retained":
                        retained = {k: event[k] for k in ("streamlines", "boundaryLayerCurves") if k in event}
                except (ValueError, AttributeError):
                    detail = line.strip() or detail
            code = process.wait()
            report.seek(0)
            try:
                final = json.load(report)
                if "residualBound" in final:
                    final["cp"] = panel_pressure(elements, final["cp"])
                final.update(retained)
            except ValueError as exc:
                if code == 0:
                    detail = f"Invalid solver report: {exc}"
                final = {}
            events.finish(code, final, detail)

    threading.Thread(target=collect, daemon=True).start()
    return process


class Plots(tk.Canvas):
    def __init__(self, parent):
        super().__init__(parent, bg=PAPER, highlightthickness=0, width=610, height=460)
        self.elements, self.frame, self.history = [], {}, []
        self.geometry_plot = tk.Canvas(self, bg=PAPER, highlightthickness=0)
        self.quantity = tk.StringVar(value=PRESSURE_VIEW)
        self.selector = ttk.Combobox(self, textvariable=self.quantity, values=(PRESSURE_VIEW, BL_VIEW),
                                     state="readonly", width=22)
        self.selector.bind("<<ComboboxSelected>>", lambda _: self.draw())
        self.bind("<Configure>", lambda _: self.draw())

    def draw(self):
        self.delete("all")
        width, height = self.winfo_width(), self.winfo_height()
        if width < 100 or height < 160:
            return
        left, right = 59, width-22
        top = 32
        geom_end = int(height*.23)
        direct = "residualBound" in self.frame
        cp_end = height-8 if direct else int(height*.70)

        def title(y, text):
            self.create_text(18, y, text=text, fill=INK, font=("TkDefaultFont", 11, "bold"), anchor="w")

        title(16, "Airfoil and flow")
        plot_w, plot_h = right-left, geom_end-top
        geom = self.geometry_plot
        geom.place(x=left, y=top, width=plot_w, height=plot_h)
        geom.delete("all")
        points = [p for _, poly in self.elements for p in poly]
        if points:
            xmin, xmax = min(p[0] for p in points), max(p[0] for p in points)
            ymin, ymax = min(p[1] for p in points), max(p[1] for p in points)
            scale = min(.82*plot_w / max(xmax-xmin, 1e-12), .73*plot_h / max(ymax-ymin, 1e-12))

            def coordinates(poly):
                return [v for x, y in poly for v in (plot_w/2 + (x-(xmin+xmax)/2)*scale,
                                                     plot_h/2 - (y-(ymin+ymax)/2)*scale)]

            for line in self.frame.get("streamlines", []):
                if len(line) >= 2:
                    geom.create_line(coordinates(line), fill="#86b8b7")
            for _, poly in self.elements:
                geom.create_polygon(coordinates(poly), fill="#d7eae7", outline=UPPER, width=1.6)
        bl = self.frame.get("boundaryLayerCurves")
        if bl:
            self.selector.place(x=16, y=geom_end+19, anchor="w")
        else:
            self.selector.place_forget()
            title(geom_end+19, PRESSURE_VIEW)
        self.create_text(right, geom_end+19, text="Upper  ━", fill=UPPER, anchor="e")
        self.create_text(right-74, geom_end+19, text="Lower  ━", fill=LOWER, anchor="e")
        colors = {"upper": UPPER, "lower": LOWER, "base": MUTED}
        if bl and self.quantity.get() == BL_VIEW:
            curves = [(b["samples"], colors[b["surface"]]) for b in bl]
            self.axes(left, geom_end+48, right, cp_end-24, curves, "s/c from stagnation", "δ*/c", origin=True)
        else:
            groups = {}
            for p in self.frame.get("cp", []):
                if p.get("cp") is not None and p.get("x") is not None:
                    groups.setdefault((p["contour"], p["surface"]), []).append((p["x"], p["cp"]))
            curves = [(pts, colors[side]) for (_, side), pts in groups.items()]
            self.axes(left, geom_end+48, right, cp_end-24, curves, "x/c", "Cp", reverse=True)
        if direct:
            return
        title(cp_end+17, "Convergence")
        curves = []
        for index, (residual, new_stage) in enumerate(self.history):
            if new_stage or not curves:
                curves.append(([], UPPER))
            curves[-1][0].append((index+1, math.log10(max(residual, 1e-16))))
        self.axes(left, cp_end+40, right, height-31, curves, "Solve step", "Residual", log=True)

    def axes(self, left, top, right, bottom, curves, xlabel, ylabel, reverse=False, log=False, origin=False):
        points = [p for line, _ in curves for p in line]
        xmin, xmax = (min(x for x, _ in points), max(x for x, _ in points)) if points else (0, 1)
        if reverse and self.elements:
            chord_x = [x for _, poly in self.elements for x, _ in poly]
            xmin, xmax = min(xmin, min(chord_x)), max(xmax, max(chord_x))
        ymin, ymax = (min(y for _, y in points), max(y for _, y in points)) if points else ((-8, 0) if log else (-2, 1))
        if origin:
            xmin = ymin = 0
        if xmin == xmax:
            xmax = xmin+1
        if ymax-ymin < 1e-8:
            ymin, ymax = ymin-.5, ymax+.5
        if log:
            ymin, ymax = math.floor(ymin), math.ceil(ymax)
        else:
            span = ymax-ymin
            ymin, ymax = (0 if origin else ymin-span*.05), ymax+span*.05

        def px(v):
            return left + (v-xmin)/(xmax-xmin)*(right-left)

        def py(v):
            t = (v-ymin)/(ymax-ymin)
            return top+t*(bottom-top) if reverse else bottom-t*(bottom-top)

        for i in range(5):
            x = xmin+(xmax-xmin)*i/4
            if not log or i in (0, 4) or xmax-xmin >= 4:
                if log:
                    x = round(x)
                label = str(x) if log else f"{x:.3g}"
                self.create_text(px(x), bottom+4, text=label, anchor="ne" if i == 4 else "n", fill=MUTED)
        yticks = (list(range(ymin, ymax+1, max(1, math.ceil((ymax-ymin)/4)))) if log
                  else [ymin+(ymax-ymin)*i/4 for i in range(5)])
        for y in yticks:
            self.create_line(left, py(y), right, py(y), fill=GRID)
            self.create_text(left-5, py(y), text=f"1e{y}" if log else f"{y:.2g}", anchor="e", fill=MUTED)
        self.create_line(left, top, left, bottom, right, bottom, fill=MUTED)
        self.create_text(left-5, top-5, text=ylabel, anchor="se", fill=MUTED)
        self.create_text((left+right)/2, bottom+18, text=xlabel, anchor="n", fill=MUTED)
        for line, color in curves:
            coords = [v for x, y in line for v in (px(x), py(y))]
            if len(coords) >= 4:
                self.create_line(coords, fill=color, width=1.8)
            elif coords:
                x, y = coords
                self.create_oval(x-2, y-2, x+2, y+2, fill=color, outline="")


class Desktop:
    def __init__(self, root):
        self.root, self.process = root, None
        self.cancelled, self.custom = False, None
        root.title("Alula • Native airfoil solver")
        root.geometry("1030x730")
        root.minsize(810, 630)
        root.configure(bg=BG)
        style = ttk.Style(root)
        style.theme_use("clam")
        style.configure("TFrame", background=BG)
        style.configure("TLabel", background=BG, foreground=INK)
        style.configure("Run.TButton", background=UPPER, foreground=PAPER, padding=(18, 7))
        style.map("Run.TButton", background=[("active", "#085d58")])
        header = tk.Frame(root, bg=INK, padx=20, pady=8)
        header.pack(fill="x")
        tk.Label(header, text="alula", font=("TkDefaultFont", 22, "bold"), bg=INK, fg=PAPER).pack(side="left")
        tk.Label(header, text="Native airfoil solver", bg=INK, fg="#b6d4d5", padx=14).pack(side="left")
        main = ttk.Frame(root, padding=12)
        main.pack(expand=True, fill="both")
        main.columnconfigure(1, weight=1)
        main.rowconfigure(0, weight=1)
        settings = ttk.Frame(main, padding=(3, 3, 16, 0))
        settings.grid(row=0, column=0, sticky="ns")
        self.controls = []
        self.case, self.mode = tk.StringVar(value=PRESETS[0]), tk.StringVar(value="Euler + boundary layer")
        ttk.Label(settings, text="Airfoil assembly", font=("TkDefaultFont", 12, "bold")).grid(row=0, columnspan=2, sticky="w", pady=(0, 9))
        self.selector = ttk.Combobox(settings, textvariable=self.case, values=PRESETS, state="readonly", width=23)
        self.selector.grid(row=1, columnspan=2, sticky="ew")
        self.selector.bind("<<ComboboxSelected>>", lambda _: self.select())
        load = ttk.Button(settings, text="Open contours…", command=self.open_contours)
        load.grid(row=2, columnspan=2, sticky="ew", pady=(5, 18))
        ttk.Label(settings, text="Flow conditions", font=("TkDefaultFont", 12, "bold")).grid(row=3, columnspan=2, sticky="w", pady=(0, 8))
        solver = ttk.Combobox(settings, textvariable=self.mode, values=("Euler + boundary layer", "Euler", PANEL_MODE, VISCOUS_MODE), state="readonly")
        solver.grid(row=4, columnspan=2, sticky="ew", pady=(0, 6))
        solver.bind("<<ComboboxSelected>>", lambda _: self.sync_controls())
        self.controls += [(self.selector, "readonly"), (load, "normal"), (solver, "readonly")]
        self.fields, self.bl_fields, self.euler_fields = [], [], []
        self.zero_mach = tk.StringVar(value="0")
        for i, name in enumerate(FIELDS):
            row = 5+i
            if i >= 6:
                row += 1
            ttk.Label(settings, text=name).grid(row=row, column=0, sticky="w", pady=4)
            var = tk.StringVar()
            self.fields.append(var)
            entry = ttk.Entry(settings, textvariable=var, width=10, justify="right")
            entry.grid(row=row, column=1, sticky="e")
            self.controls.append((entry, "normal"))
            if i == 1:
                self.mach_entry = entry
            if 2 <= i < 6:
                self.bl_fields.append(entry)
            if i != 0:
                self.euler_fields.append(entry)
        ttk.Label(settings, text="Mesh", font=("TkDefaultFont", 12, "bold")).grid(row=11, columnspan=2, sticky="w", pady=(14, 4))
        self.action = ttk.Button(settings, text="Run", style="Run.TButton", command=self.run_or_cancel)
        self.action.grid(row=14, columnspan=2, sticky="ew", pady=(18, 4))
        self.status = tk.StringVar()
        ttk.Label(settings, textvariable=self.status, wraplength=203, justify="left").grid(row=15, columnspan=2, sticky="nw", pady=6)
        visual = ttk.Frame(main)
        visual.grid(row=0, column=1, sticky="nsew")
        self.plots = Plots(visual)
        self.plots.pack(fill="both", expand=True)
        self.loads = tk.Label(visual, bg=INK, fg=PAPER, padx=12, pady=10, anchor="w")
        self.loads.pack(fill="x")
        root.protocol("WM_DELETE_WINDOW", self.close)
        self.select()

    def select(self):
        name = self.case.get()
        try:
            self.plots.elements = preset_geometry(name) if name in PRESETS else self.custom
        except (OSError, ValueError) as exc:
            self.plots.elements = []
            self.status.set(str(exc))
            return
        defaults = DEFAULTS[PRESETS.index(name)] if name in PRESETS else DEFAULTS[0]
        for var, value in zip(self.fields, defaults):
            var.set(value)
        self.reset_plots()
        self.status.set("Ready")

    def reset_plots(self):
        self.plots.frame, self.plots.history = {}, []
        self.last_stage = None
        self.update_loads({})
        self.plots.draw()

    def open_contours(self):
        path = filedialog.askopenfilename(title="Open native contours", filetypes=(("Contour files", "*.txt"), ("All files", "*")))
        if not path:
            return
        try:
            self.custom = read_contours(Path(path))
        except (OSError, ValueError) as exc:
            self.status.set(str(exc))
            return
        label = "File · " + Path(path).name
        self.selector.configure(values=(*PRESETS, label))
        self.case.set(label)
        self.select()

    def update_loads(self, frame):
        def number(k):
            value = frame.get(k)
            return f"{value:.5f}" if isinstance(value, (int, float)) else "—"
        self.loads.config(text="    ".join(f"{label}  {number(k)}" for label, k in (
            ("CL", "cl"), ("Cm", "cm"), ("CD pressure", "pressureDrag"), ("CD wake", "viscousWakeDrag"))))

    def sync_controls(self):
        mode = self.mode.get()
        self.mach_entry.configure(textvariable=self.zero_mach if mode in (PANEL_MODE, VISCOUS_MODE) else self.fields[1])
        for widget, state in self.controls:
            disabled = (self.process is not None or (mode == "Euler" and widget in self.bl_fields)
                        or (mode == PANEL_MODE and widget in self.euler_fields)
                        or (mode == VISCOUS_MODE and widget in self.euler_fields and widget not in self.bl_fields))
            widget.configure(state="disabled" if disabled else state)

    def run_or_cancel(self):
        if self.process is not None:
            self.cancelled = True
            self.process.kill()
            self.status.set("Cancelling…")
            self.action.configure(state="disabled")
            return
        self.events = Updates()
        try:
            if not self.plots.elements:
                raise ValueError("Select an airfoil or open a contour file.")
            args = solver_args(self.case.get(), self.mode.get(), [v.get() for v in self.fields], self.plots.elements)
            self.process = launch(args, self.plots.elements, self.events)
        except (OSError, ValueError) as exc:
            self.status.set(str(exc))
            return
        self.cancelled = False
        self.reset_plots()
        self.status.set("Solving panel system…" if self.mode.get() in (PANEL_MODE, VISCOUS_MODE) else "Building mesh…")
        self.action.configure(text="Cancel")
        self.sync_controls()
        self.root.after(20, self.poll)

    def poll(self):
        samples, latest, finished = self.events.take()
        for name, mach, step, residual in samples:
            stage = (name, mach)
            if residual is not None:
                self.plots.history.append((residual, stage != self.last_stage or step == 0))
            self.last_stage = stage

        if finished is not None:
            code, final, detail = finished
            if code == 0 and "residualBound" in final:
                status = f"Panel solve complete\nResidual bound {final['residualBound']:.3g}"
            elif code == 0 and final.get("converged"):
                status = "Converged"
            elif self.cancelled:
                status = "Cancelled"
                if latest or self.plots.frame:
                    status += " · last iterate"
            elif final.get("model") == "panel-bl":
                status = final["reason"].capitalize()
            elif final:
                status = "Did not converge" if final.get("gridStatus") == 0 else "Mesh generation failed"
            else:
                status = detail or f"Solver exited ({code})"
            if final:
                latest = final
                coupled = final.get("coupled", {})
                norms = [v for v in (final.get("residual"), coupled.get("edgeResidual"),
                                     coupled.get("blResidual"), coupled.get("displacementResidual")) if v is not None]
                if final.get("equationResidual") is not None:
                    status += f"\nEquation residual {final['equationResidual']:.3g}"
                    if final.get("wakeResidual") is not None:
                        status += f"\nWake residual {final['wakeResidual']:.3g}"
                elif norms:
                    status += f"\nResidual {max(norms):.3g}"
                actual = final.get("actualMach", final.get("mach"))
                if actual is not None:
                    status += f" · M {actual:.4g}"
            self.status.set(status)
            self.process = None
            self.action.configure(text="Run", state="normal")
            self.sync_controls()
        elif latest is not None and not self.cancelled:
            residual = latest["residual"]
            incidence = f" · α {latest['actualAlpha']:.3g}°" if "actualAlpha" in latest else ""
            if latest.get("wakeResidual") is not None and latest.get("equationResidual") is not None:
                incidence += f"\nEquation residual {latest['equationResidual']:.3g}\nWake residual {latest['wakeResidual']:.3g}"
                residual = None
            self.status.set(f"{latest['stage']} · iteration {latest['iteration']}\nM {latest['mach']:.4g}{incidence}" +
                            (f" · residual {residual:.3g}" if residual is not None else ""))
        if latest is not None:
            self.plots.frame = latest
            self.update_loads(latest)
            self.plots.draw()
        if self.process is not None:
            self.root.after(30, self.poll)

    def close(self):
        if self.process is not None:
            self.process.kill()
            self.process.wait()
        self.root.destroy()


if __name__ == "__main__":
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        raise SystemExit(f"Could not open the desktop window: {exc}") from None
    Desktop(root)
    root.mainloop()
