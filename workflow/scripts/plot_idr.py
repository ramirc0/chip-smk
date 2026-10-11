"""Redraw the `idr --plot` panels from the IDR output.

Reads every merged peak from the narrowPeak that idr writes with
`--soft-idr-threshold`, so IDR itself never plots. The four panels match
idr 2.0.4: rank vs rank and score vs score, colored by the IDR threshold, then
-log10 IDR against each sample's rank in 20 bins. idr writes -log10 IDR capped
at 5, so the top bins clip there.
"""

import argparse
import math

import matplotlib.pyplot as plt
import numpy as np
import poikilos as pk
import polars as pl


def build_parser():
    """Build the plot_idr.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_idr.py flags.
    """
    p = argparse.ArgumentParser(description="Plot IDR results.")
    p.add_argument(
        "-i", "--input", required=True, help="Gzipped idr narrowPeak output."
    )
    p.add_argument("--threshold", type=float, default=0.05, help="IDR threshold.")
    p.add_argument("--title", default="", help="Figure title.")
    p.add_argument(
        "-o", "--output", required=True, help="SVG to write (PNG beside it)."
    )
    return p


def main(argv=None):
    """Draw the four IDR panels, or a note when IDR did not run.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    peaks = _read(args.input)

    pk.use("plain")
    if peaks.is_empty():
        fig, ax = plt.subplots(figsize=(4, 1.5))
        ax.text(0.5, 0.5, "Too few peaks for IDR", ha="center", va="center")
        ax.set_axis_off()
    else:
        fig, axes = plt.subplots(2, 2, figsize=(8, 8))
        _draw(axes, peaks, args.threshold)
    fig.suptitle(args.title)
    pk.save_figure(fig, args.output)
    plt.close(fig)


def _read(path):
    """Read the IDR columns and rank each sample's scores.

    Column 12 is the global -log10 IDR. Columns 15 and 19 are each sample's
    score that IDR ranked on.

    Parameters
    ----------
    path : str
        Gzipped idr narrowPeak output. Empty when IDR had too few peaks.

    Returns
    -------
    polars.DataFrame
        `neg_log10_idr`, `score_*`, and `rank_*` in (0, 1], highest score last.
    """
    try:
        t = pl.read_csv(path, separator="\t", has_header=False)
    except pl.exceptions.NoDataError:
        return pl.DataFrame()
    t = t.select(
        pl.nth(11).alias("neg_log10_idr"),
        pl.nth(14).alias("score_1"),
        pl.nth(18).alias("score_2"),
    )
    # idr breaks ties at random; a fixed seed keeps the plot reproducible.
    return t.with_columns(
        (pl.col(s).rank("random", seed=0) / pl.len()).alias(f"rank_{s[-1]}")
        for s in ("score_1", "score_2")
    )


def _draw(axes, peaks, threshold):
    """Draw the four panels.

    Parameters
    ----------
    axes : numpy.ndarray of matplotlib.axes.Axes
        2 x 2 axes.
    peaks : polars.DataFrame
        Output of `_read`.
    threshold : float
        IDR threshold. Peaks above it are drawn in red.
    """
    colors = pk.PALETTES["plain"].colors
    passed = peaks["neg_log10_idr"].to_numpy() >= -math.log10(threshold)
    color = np.where(passed, pk.PALETTES["plain"].foreground, colors["red"])
    title = f"red: IDR > {threshold:g}"
    dots = {"s": 4, "alpha": 0.25, "linewidths": 0, "rasterized": True}

    ax = axes[0, 0]
    ax.scatter(peaks["rank_1"], peaks["rank_2"], c=color, **dots)
    ax.set(xlabel="Sample 1 rank", ylabel="Sample 2 rank", title=f"Ranks ({title})")
    ax = axes[0, 1]
    ax.scatter(
        np.log10(peaks["score_1"] + 1), np.log10(peaks["score_2"] + 1), c=color, **dots
    )
    ax.set(
        xlabel="Sample 1 log10 score",
        ylabel="Sample 2 log10 score",
        title=f"Scores ({title})",
    )
    for ax in axes[0]:
        pk.despine(ax)

    for i, ax in enumerate(axes[1], start=1):
        rank = peaks[f"rank_{i}"].to_numpy()
        idr = peaks["neg_log10_idr"].to_numpy()
        bins = np.minimum((rank * 20).astype(int), 19)
        ax.scatter(
            rank,
            idr,
            color=colors["gray"],
            s=4,
            alpha=0.1,
            linewidths=0,
            rasterized=True,
        )
        ax.boxplot(
            [idr[bins == b] for b in range(20)],
            positions=(np.arange(20) + 0.5) / 20,
            widths=0.03,
            showfliers=False,
            manage_ticks=False,
        )
        ax.set(
            xlabel=f"Sample {i} rank",
            ylabel="-log10 IDR",
            title=f"Sample {i} rank vs IDR",
        )
        pk.despine(ax)


if __name__ == "__main__":
    main()
