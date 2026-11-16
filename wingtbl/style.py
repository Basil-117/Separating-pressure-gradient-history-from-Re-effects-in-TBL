"""Shared plotting style.

The four wing cases form an ordered Reynolds-number sequence, so they get a
single-hue ordinal ramp (light -> dark = low Re -> high Re) rather than four
unrelated categorical hues: the colour then carries the trend that the figures
are about. Validated with the data-viz palette validator in ``--ordinal`` mode
against the light surface (monotone lightness, adjacent dL >= 0.06, light end
2.06:1 vs surface, hue spread 4 deg).

W1 is drawn dashed throughout. The paper shows it is not a well-behaved TBL --
too low a Reynolds number to sustain a proper turbulent state under this
pressure gradient -- so it should never be read as just "the lowest curve".
That also gives identity a second channel besides colour.
"""

from __future__ import annotations

import matplotlib as mpl

SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e4e3df"

# Ordinal blue ramp, steps 250 / 400 / 550 / 700.
CASE_COLOR = {
    "W1": "#86b6ef",
    "W2": "#3987e5",
    "W4": "#1c5cab",
    "W10": "#0d366b",
}

CASE_LABEL = {
    "W1": r"$Re_c=100\,$k",
    "W2": r"$Re_c=200\,$k",
    "W4": r"$Re_c=400\,$k",
    "W10": r"$Re_c=1\,$M",
}

# W1 is not a well-behaved TBL; see module docstring.
CASE_DASH = {"W1": (3, 2), "W2": None, "W4": None, "W10": None}

CASE_ORDER = ["W1", "W2", "W4", "W10"]

ZPG_COLOR = "#52514e"


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "axes.edgecolor": TEXT_SECONDARY,
            "axes.labelcolor": TEXT_PRIMARY,
            "axes.titlecolor": TEXT_PRIMARY,
            "axes.linewidth": 0.8,
            "axes.grid": True,
            "axes.axisbelow": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "xtick.color": TEXT_SECONDARY,
            "ytick.color": TEXT_SECONDARY,
            "xtick.labelcolor": TEXT_SECONDARY,
            "ytick.labelcolor": TEXT_SECONDARY,
            "text.color": TEXT_PRIMARY,
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.titleweight": "semibold",
            "legend.frameon": False,
            "legend.fontsize": 8.5,
            "lines.linewidth": 2.0,
            "lines.markersize": 4,
            "figure.dpi": 130,
            "savefig.dpi": 200,
            "savefig.bbox": "tight",
        }
    )


def style_axes(ax) -> None:
    """Recede the frame: keep left and bottom spines only."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(TEXT_SECONDARY)
