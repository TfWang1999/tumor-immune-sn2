"""Approved three-panel Figure 3 in manuscript dimensionless coordinates.

Self-contained; only saved CSV/JSON and algebraic theory expressions.
alpha=lambdaT/r, sigma=S/r**2, C=C_tau0/rho**3, C0=C0_tau0/rho.
No roots, model trajectories, continuation or coefficient fits are computed.
Run: python plot_figure_03.py
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


def figure_03():
    rows = read_csv(FOLLOW / "sn2_weak_recruitment.csv")
    derived = []
    for row in rows:
        sigma = number(row, "S_per_day2") / GROWTH**2
        alpha = number(row, "lambdaT_SN2_per_day") / GROWTH
        leading = number(row, "leading_prediction_per_day") / GROWTH
        derived.append({
            "model": row["model"], "sigma": sigma, "alpha": alpha,
            "residual": alpha - leading,
            "C": number(row, "C_exact") / RHO**3,
            "C0": number(row, "C0_exact") / RHO,
        })
    differences = []
    by_key = {(row["model"], row["sigma"]): row for row in derived}
    for sigma in sorted({row["sigma"] for row in derived}):
        m3, m4 = by_key[("M3", sigma)], by_key[("M4", sigma)]
        differences.append({"sigma": sigma, "difference": m4["alpha"] - m3["alpha"],
                            "limit": m4["C0"] - m3["C0"]})

    fig = plt.figure(figsize=(180 * MM, 95 * MM))
    a = fig.add_axes((0.09, 0.33, 0.23, 0.50))
    b = fig.add_axes((0.40, 0.33, 0.23, 0.50))
    c = fig.add_axes((0.71, 0.33, 0.23, 0.50))
    ss_dense = np.geomspace(1e-5 / GROWTH**2, 1.28e-3 / GROWTH**2, 120)
    for model, marker in (("M3", "o"), ("M4", "s")):
        table = sorted((r for r in derived if r["model"] == model), key=lambda r: r["sigma"])
        ss = np.array([r["sigma"] for r in table])
        ll = np.array([r["alpha"] for r in table])
        residual = np.array([r["residual"] for r in table]) * 10
        color = COLORS[model]
        # Theoretical lines evaluate existing coefficients algebraically.
        # Lines joining saved points in b/c only guide the eye. No new roots
        # or continuation points are computed.
        a.plot(ss_dense, table[0]["C"] / ss_dense + table[0]["C0"],
               color=color, ls=(0, (4, 2.5)), lw=1.0)
        a.plot(ss, ll, marker=marker, ls="none", color=color, mfc="white",
               ms=5 if model == "M3" else 4.5, mew=1.1, zorder=5)
        b.plot(ss, residual, marker=marker, color=color, mfc="white", ms=4.5,
               mew=1.0, lw=1.0, zorder=4)
        b.axhline(table[0]["C0"] * 10, color=color, ls=(0, (4, 2.5)), lw=1.0)
        b.text(1.12e-5 / GROWTH**2, table[0]["C0"] * 10 + (0.021 if model == "M3" else -0.049),
               r"$10C_{0," + model[-1] + r"}$", color=color, fontsize=9)

    C = derived[0]["C"]
    a.plot(ss_dense, C / ss_dense, color=GREY, ls=(0, (1.2, 2)), lw=1.4, zorder=2)
    a.set_yscale("log")
    a.set(ylim=(1, 300), ylabel=r"SN2 threshold $\alpha_{\mathrm{SN2}}$")
    a.yaxis.set_major_locator(LogLocator(base=10, numticks=4))
    a.yaxis.set_minor_locator(LogLocator(base=10, subs=(2, 5), numticks=12))
    a.yaxis.set_minor_formatter(NullFormatter())
    a.text(0.09, 0.12, r"$\alpha_{\mathrm{SN2}}\sim C/\sigma$", transform=a.transAxes, fontsize=9.5)
    b.set(ylim=(2.07, 2.49), ylabel=r"$10(\alpha_{\mathrm{SN2}}-C/\sigma)$")
    b.yaxis.set_major_locator(FixedLocator((2.1, 2.2, 2.3, 2.4)))

    ss = np.array([r["sigma"] for r in differences])
    delta = np.array([r["difference"] for r in differences]) * 100
    limit = differences[0]["limit"] * 100
    c.plot(ss, delta, color=INK, marker="D", mfc="white", mew=1.0, ms=4.5, lw=1.0)
    c.axhline(limit, color=INK, ls=(0, (4, 2.5)), lw=1.0)
    c.text(1.15e-5 / GROWTH**2, limit + 0.022, r"$10^2\delta\gamma/(2\beta)$", fontsize=9.5)
    c.set(ylim=(2.50, 2.765), ylabel=r"$10^2\Delta\alpha$")
    c.yaxis.set_major_locator(FixedLocator((2.5, 2.55, 2.6, 2.65, 2.7)))

    for ax, letter, title in ((a, "a", "Weak-recruitment scaling"),
                             (b, "b", "Constant correction"),
                             (c, "c", "M4 minus M3")):
        style_axis(ax, log_x=True)
        ax.set_xlim(7.5e-6 / GROWTH**2, 1.8e-3 / GROWTH**2)
        ax.set_xlabel(r"Recruitment $\sigma$")
        panel_title(ax, letter, title)
    handles = [
        Line2D([], [], color=TEAL, ls="none", marker="o", mfc="white", mew=1.1, ms=5, label="M3: linear DC"),
        Line2D([], [], color=PURPLE, ls="none", marker="s", mfc="white", mew=1.1, ms=4.5, label="M4: saturated DC"),
        Line2D([], [], color=INK, ls=(0, (4, 2.5)), lw=1, label="Asymptotic expressions"),
        Line2D([], [], color=GREY, ls=(0, (1.2, 2)), lw=1.4, label="Leading term only (a)"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.52, 0.022),
               ncol=2, frameon=False, columnspacing=2.4, handlelength=2.2, labelspacing=0.6)
    save_figure(fig, "Figure_03")


if __name__ == "__main__":
    figure_03()
