"""Approved four-panel Figure 2 in manuscript dimensionless coordinates.

Self-contained; read saved results only. alpha=lambdaT/r, gamma=dD/r,
sigma=S/r**2. State x and stability flags are unchanged. No numerical
model, continuation or other figure script is imported.
Run: python plot_figure_02.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, LogLocator, NullFormatter
import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
ROOT = HERE
MECH = ROOT / "data" / "equilibria"
FOLLOW = ROOT / "data" / "sn2"
CONFIG = ROOT / "data" / "parameters.json"
PARAMETERS = json.loads(CONFIG.read_text(encoding="utf-8"))["analysis"]
GROWTH = float(PARAMETERS["dimensional"]["r"])
RHO = GROWTH / float(PARAMETERS["r0_per_day"])

MM = 1 / 25.4
TEAL = "#168383"
PURPLE = "#74429A"
INK = "#30343B"
GREY = "#858B91"
LIGHT = "#DDE1E5"
COLORS = {"M1": TEAL, "M2": PURPLE, "M3": TEAL, "M4": PURPLE}

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 10,
        "axes.titlesize": 10.5,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "axes.linewidth": 0.75,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "lines.linewidth": 1.35,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "svg.hashsalt": "tumour-immune-dimensionless",
        "savefig.facecolor": "white",
        "savefig.bbox": None,
    }
)


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def number(row, name):
    return float(row[name])


def is_stable(row):
    return row["stable"] == "True"


def style_axis(ax, *, log_x=False):
    ax.set_axisbelow(True)
    ax.grid(axis="y", which="major", color=LIGHT, linewidth=0.5, alpha=0.65)
    ax.tick_params(direction="out", length=3, width=0.65, pad=3)
    ax.spines[["top", "right"]].set_visible(False)
    if log_x:
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(LogLocator(base=10, numticks=4))
        ax.xaxis.set_minor_locator(LogLocator(base=10, subs=(2, 5), numticks=12))
        ax.xaxis.set_minor_formatter(NullFormatter())


def panel_title(ax, letter, title):
    ax.set_title(title, loc="left", pad=10)
    ax.text(
        -0.16 if ax.get_position().width < 0.5 else -0.075,
        1.05,
        letter,
        transform=ax.transAxes,
        fontweight="bold",
        fontsize=12,
        va="bottom",
    )


def draw_stability_paths(ax, rows, x_column, y_column, color, *, x_factor=1.0):
    """Mask opposite stability with NaN, never reconnect across removed rows.

    Model/kind selection has already identified a single equilibrium curve.
    Its order is a monotonic state-grid order, NOT a sorted parameter order.
    A fold therefore remains a turning point. Independent model/kind curves
    are separate calls and can never be joined together.
    """
    horizontal = np.array([number(r, x_column) * x_factor for r in rows])
    vertical = np.array([number(r, y_column) for r in rows])
    flags = np.array([is_stable(r) for r in rows])
    for stable, line_style in ((True, "-"), (False, (0, (4, 2.5)))):
        ordinate = np.where(flags == stable, vertical, np.nan)
        ax.plot(horizontal, ordinate, color=color, ls=line_style, zorder=2)


def save_figure(fig, stem):
    OUT.mkdir(exist_ok=True)
    for extension in ("pdf", "svg", "png"):
        options = {"dpi": 300} if extension == "png" else {}
        fig.savefig(OUT / f"{stem}.{extension}", **options)
    plt.close(fig)
    print(f"Saved {stem}: PDF, SVG, PNG", flush=True)


def main():
    OUT.mkdir(exist_ok=True)
    rows = read_csv(MECH / "global_equilibrium_branches.csv")
    summary = {r["model"]: r for r in read_csv(MECH / "ablation_summary.csv")}
    landmarks = read_csv(MECH / "ablation_landmarks.csv")
    folds = {r["model"]: r for r in read_csv(MECH / "upper_fold_checks.csv")}

    fig = plt.figure(figsize=(180 * MM, 165 * MM))
    a = fig.add_axes((0.105, 0.630, 0.365, 0.295))
    b = fig.add_axes((0.605, 0.630, 0.365, 0.295))
    c = fig.add_axes((0.105, 0.185, 0.365, 0.295))
    d = fig.add_axes((0.605, 0.185, 0.365, 0.295))

    def internal_curves(ax, models):
        for model in models:
            curve = [r for r in rows if r["model"] == model and r["kind"] == "internal"]
            # Same curve/order in global and local panels. View limits crop
            # the line; no out-of-window points are deleted and reconnected.
            draw_stability_paths(ax, curve, "lambdaT_per_day", "x", COLORS[model], x_factor=1/GROWTH)

    def boundary_high(ax):
        curve = [r for r in rows if r["model"] == "M1" and r["kind"] == "boundary_tumour"]
        draw_stability_paths(ax, curve, "lambdaT_per_day", "x", INK, x_factor=1/GROWTH)

    for ax, models, letter, title in (
        (a, ("M1", "M2"), "a", "Global: without recruitment"),
        (b, ("M3", "M4"), "b", "Global: with recruitment"),
    ):
        internal_curves(ax, models)
        for model in models:
            curve = [r for r in rows if r["model"] == model and r["kind"] == "internal"]
            saved_sn = next((r for r in landmarks if r["model"] == model
                             and r["kind"] == "equilibrium" and r["type"] == "LP"), None)
            sn_x = number(saved_sn, "x_min") if saved_sn else number(
                min((r for r in curve if number(r, "x") < 0.9),
                    key=lambda r: number(r, "lambdaT_per_day")), "x")
            ax.plot(number(summary[model], "lambdaT_SN_per_day") / GROWTH, sn_x, "D", ms=4.5,
                    mfc="white", mec=COLORS[model], mew=1, zorder=5)
            hb = next(r for r in landmarks if r["model"] == model
                      and r["kind"] == "equilibrium" and r["type"] == "HB")
            ax.plot(number(hb, "lambdaT_per_day") / GROWTH, number(hb, "x_min"), "o", ms=4.5,
                    mfc="white", mec=COLORS[model], mew=1, zorder=5)

        origin = [r for r in rows if r["model"] == models[0] and r["kind"] == "boundary_origin"]
        draw_stability_paths(ax, origin, "lambdaT_per_day", "x", GREY, x_factor=1/GROWTH)
        if ax is a:
            boundary_high(ax)
            ax.plot(0.84 / GROWTH, 1, "s", color=INK, mfc="white", ms=4.7, mew=1)
            ax.text(0.91 / GROWTH, 0.94, r"$E_K$", fontsize=10)
            ax.text(0.52 / GROWTH, 0.79, "M1: linear DC\nM2: saturated DC", fontsize=9,
                    va="top", linespacing=1.35)
        else:
            for model in models:
                fold = folds[model]
                ax.plot(number(fold, "lambdaT_upper_SN_per_day") / GROWTH, number(fold, "x"), "D",
                        ms=4.5, mfc="white", mec=COLORS[model], mew=1, zorder=5)
            ax.annotate(r"$E_H$", xy=(0.64 / GROWTH, 0.997), xytext=(0.69 / GROWTH, 0.86),
                        arrowprops={"arrowstyle": "-", "lw": 0.65, "color": INK}, fontsize=10)
            ax.annotate("SN2", xy=(1.203 / GROWTH, 0.99536), xytext=(1.0 / GROWTH, 0.77),
                        arrowprops={"arrowstyle": "-", "lw": 0.65, "color": INK}, fontsize=9.5)
            ax.text(0.48 / GROWTH, 0.66, "M3: linear DC\nM4: saturated DC", fontsize=9,
                    va="top", linespacing=1.35)
        ax.annotate("SN1", xy=(0.17 / GROWTH, 0.5), xytext=(0.32 / GROWTH, 0.60),
                    arrowprops={"arrowstyle": "-", "lw": 0.65, "color": GREY}, fontsize=9.5)
        ax.annotate("Hopf", xy=(0.185 / GROWTH, 0.251), xytext=(0.36 / GROWTH, 0.29),
                    arrowprops={"arrowstyle": "-", "lw": 0.65, "color": GREY}, fontsize=9.5)
        ax.text(1.10 / GROWTH, 0.055, r"$O$", color=GREY, fontsize=10)
        ax.set(xlim=(0.10 / GROWTH, 1.25 / GROWTH), ylim=(-0.025, 1.065),
               xlabel=r"T-cell expansion $\alpha$",
               ylabel=r"Tumour fraction $x=C/K$")
        ax.xaxis.set_major_locator(FixedLocator((0.5, 1, 2, 3)))
        ax.yaxis.set_major_locator(FixedLocator((0, 0.25, 0.5, 0.75, 1)))
        style_axis(ax)
        panel_title(ax, letter, title)

    internal_curves(c, ("M1", "M2"))
    boundary_high(c)
    # Use an actual saved boundary point for the shared M1/M2 state marker.
    ek = min((r for r in rows if r["model"] == "M1" and r["kind"] == "boundary_tumour"),
             key=lambda r: abs(number(r, "lambdaT_per_day") - 1.1))
    c.plot(number(ek, "lambdaT_per_day") / GROWTH, number(ek, "x"), "s", ms=4.7,
           mfc="white", mec=INK, mew=1, zorder=5)
    c.text(0.04, 0.86, r"Stable $E_K$ (M1/M2)", transform=c.transAxes, fontsize=9.5)
    c.text(0.05, 0.55, "Upper state persists", transform=c.transAxes, fontsize=9.5)
    c.text(0.50, 0.07, "Saddles", transform=c.transAxes, fontsize=9)

    internal_curves(d, ("M3", "M4"))
    for model in ("M3", "M4"):
        fold = folds[model]
        d.plot(number(fold, "lambdaT_upper_SN_per_day") / GROWTH, number(fold, "x"), "D",
               ms=4.5, mfc="white", mec=COLORS[model], mew=1, zorder=5)
    d.text(0.04, 0.83, r"Stable interior $E_H$", transform=d.transAxes, fontsize=9.5)
    d.annotate("SN2", xy=(1.203 / GROWTH, 0.99536), xytext=(1.205 / GROWTH, 0.9973),
               arrowprops={"arrowstyle": "-", "lw": 0.65, "color": INK}, fontsize=9.5)
    d.text(0.37, 0.24, "Saddles", transform=d.transAxes, fontsize=9)

    for ax, letter, title in ((c, "c", "Local: without recruitment"),
                              (d, "d", "Local: with recruitment")):
        # Identical absolute coordinates: control and recruited branches can
        # be compared directly, without assigning a fictitious control SN2.
        for model in ("M3", "M4"):
            ax.axvline(number(folds[model], "lambdaT_upper_SN_per_day") / GROWTH,
                       color=LIGHT, lw=0.65, ls=(0, (2, 3)), zorder=1)
        ax.set(xlim=(0.90 / GROWTH, 1.25 / GROWTH), ylim=(0.9855, 1.0012),
               xlabel=r"T-cell expansion $\alpha$",
               ylabel=r"Tumour fraction $x=C/K$")
        ax.xaxis.set_major_locator(FixedLocator((2.4, 2.6, 2.8, 3.0, 3.2)))
        ax.yaxis.set_major_locator(FixedLocator((0.988, 0.992, 0.996, 1.000)))
        ax.ticklabel_format(axis="y", useOffset=False, style="plain")
        style_axis(ax)
        panel_title(ax, letter, title)

    handles = [
        Line2D([], [], color=TEAL, label="Linear DC"),
        Line2D([], [], color=PURPLE, label="Saturated DC"),
        Line2D([], [], color=INK, ls="-", label="Stable equilibrium"),
        Line2D([], [], color=INK, ls=(0, (4, 2.5)), label="Unstable equilibrium"),
        Line2D([], [], color=INK, marker="D", ls="none", mfc="white", ms=4.5, label="Saddle-node"),
        Line2D([], [], color=INK, marker="o", ls="none", mfc="white", ms=4.5, label="Hopf"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.54, 0.017),
               ncol=3, frameon=False, columnspacing=1.6, handlelength=2.2, labelspacing=0.55)
    for extension in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"Figure_02.{extension}", **({"dpi": 300} if extension == "png" else {}))
    plt.close(fig)
    print("Saved revised four-panel Figure 2: PDF, SVG, PNG", flush=True)


if __name__ == "__main__":
    main()
