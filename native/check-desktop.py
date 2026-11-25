#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Exercise the installed desktop viewer and C binaries on an X display."""
import time
import tkinter as tk

from desktop import BL_VIEW, Desktop, PANEL_MODE, VISCOUS_MODE, PRESETS, panel_pressure, preset_geometry


def wait_for(root, predicate):
    deadline = time.monotonic() + 90
    while not predicate():
        root.update()
        if time.monotonic() > deadline:
            raise AssertionError("Timed out waiting for the native desktop solver")
        time.sleep(.01)


def main():
    root = tk.Tk()
    app = Desktop(root)
    try:
        app.case.set("NLR 7301 + flap")
        app.select()
        app.mode.set(PANEL_MODE)
        app.sync_controls()
        assert app.mach_entry.get() == "0" and float(app.fields[1].get()) == .185
        app.fields[1].set("")  # Panel flow only requires an angle of attack.
        app.run_or_cancel()
        assert app.process is not None
        wait_for(root, lambda: app.process is None)
        report = app.plots.frame
        assert app.status.get().startswith("Panel solve complete")
        assert report["residualBound"] < 1e-9
        assert len(report["cp"]) == 704
        assert sum(p["surface"] == "base" for p in report["cp"]) == 64
        assert not app.plots.history

        # Move the closing point onto the wetted surface. Each Cp still
        # belongs to the same physical panel and the plot order is continuous.
        shifted, values, reordered = [], [], []
        for kind, points in preset_geometry("NLR 7301 + flap"):
            n = len(points)-1
            _, upper, lower = kind.split()
            shift = 47
            shifted.append((f"F {(int(upper)-shift)%n} {(int(lower)-shift)%n}",
                            points[shift:n]+points[:shift+1]))
            samples = list(range(len(values), len(values)+n))
            values.extend(samples)
            reordered.extend(samples[shift:]+samples[:shift])
        assert panel_pressure(app.plots.elements, values) == panel_pressure(shifted, reordered)

        app.case.set(PRESETS[0])
        app.select()
        app.mode.set(VISCOUS_MODE)
        app.sync_controls()
        app.fields[1].set("")  # Incompressible viscous flow fixes Mach to zero.
        app.run_or_cancel()
        assert app.process is not None
        wait_for(root, lambda: app.plots.frame.get("cp") or app.process is None)
        assert app.plots.history and app.plots.frame["boundaryLayerCurves"] and app.plots.frame["streamlines"]
        wait_for(root, lambda: app.process is None)
        assert app.status.get().startswith("Converged")
        assert app.plots.frame["model"] == "panel-bl" and app.plots.frame["wakeResidual"] <= 1e-6
        assert app.plots.frame["equationResidual"] <= 1e-8 and app.plots.frame["viscousWakeDrag"] > 0

        app.mode.set("Euler + boundary layer")
        app.case.set("RAE 2822")
        app.select()
        app.sync_controls()
        assert float(app.mach_entry.get()) == .74
        app.fields[6].set(64)
        app.fields[7].set(24)
        app.run_or_cancel()
        child = app.process
        assert child is not None
        wait_for(root, lambda: app.plots.frame.get("cp") or app.process is None)
        assert app.process is child and app.plots.frame["streamlines"] and app.plots.history
        assert all(app.plots.frame[k] is not None for k in ("cl", "cm", "pressureDrag"))
        app.run_or_cancel()
        wait_for(root, lambda: app.process is None)
        assert app.status.get().startswith("Cancelled") and child.poll() is not None

        app.case.set(PRESETS[0])
        app.select()
        app.run_or_cancel()
        wait_for(root, lambda: app.process is None)
        assert app.status.get().startswith("Converged")
        assert app.plots.frame["viscousWakeDrag"] is not None
        assert len(app.plots.frame["boundaryLayerCurves"]) == 2
        app.plots.quantity.set(BL_VIEW)
        app.plots.selector.event_generate("<<ComboboxSelected>>")
        root.update()
        labels = [app.plots.itemcget(item, "text") for item in app.plots.find_all()
                  if app.plots.type(item) == "text"]
        assert "δ*/c" in labels and "s/c from stagnation" in labels
    finally:
        app.close()
    print("Desktop panel, viscous wake, live Euler progress, cancellation and coupled restart passed.")


if __name__ == "__main__":
    main()
