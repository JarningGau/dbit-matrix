"""Nature-style matplotlib theme (companion of ``theme_nature.R``).

Usage: copy this file next to the plotting script, or add this folder to the
import path, then apply the theme *after* any seaborn theme::

    import sys
    sys.path.insert(0, "<skill>/scripts")   # skip if the file is next to the script
    import theme_nature as tn

    tn.apply_theme()
    fig, ax = plt.subplots(figsize=(tn.NATURE_W1_IN, tn.NATURE_W1_IN * 0.7))
    ...
    tn.style_axes(ax)
    tn.save_nature_pdf(fig, "figure.pdf")

Text is 5-7 pt Arial/Helvetica, PDF text stays selectable (fonttype 42), the
background is white, there is no grid, spines and ticks are black 0.5 pt, and
the top/right spines are hidden. ``savefig.dpi`` is 600 (minimum 300).
"""
import matplotlib as mpl

## Nature Research figure widths (inches)
NATURE_W1_IN = 88 / 25.4  # 1 column
NATURE_W2_IN = 180 / 25.4  # 2 columns

RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 6,
    "axes.titlesize": 6,
    "axes.labelsize": 6,
    "xtick.labelsize": 5,
    "ytick.labelsize": 5,
    "legend.fontsize": 5,
    "figure.titlesize": 6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 600,
    "axes.linewidth": 0.5,
    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "axes.facecolor": "white",
    "axes.grid": False,
    "figure.facecolor": "white",
    "text.color": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2,
    "ytick.major.size": 2,
    "legend.frameon": False,
}


def apply_theme():
    """Install the Nature rcParams."""
    mpl.rcParams.update(RC)


def style_axes(ax):
    """Open top/right spines, black 0.5 pt spines and ticks, no grid."""
    ax.set_axisbelow(False)
    ax.grid(False)
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color("black")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(colors="black", width=0.5, length=2, pad=1, labelsize=5)


def save_nature_pdf(fig, path):
    """Save one vector PDF at the figure's drawn size (no tight bbox)."""
    fig.savefig(path)
    return path
