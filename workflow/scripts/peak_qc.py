"""Report peak counts and region size statistics, as ENCODE's peak QC does.

ENCODE summarizes widths with pandas ``describe()``: linear-interpolated
quartiles, minimum, maximum and mean.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import poikilos as pk
import polars as pl

STATS = ["min_size", "25_pct", "50_pct", "75_pct", "max_size", "mean"]


def build_parser():
    """Build the peak_qc.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the peak_qc.py flags.
    """
    p = argparse.ArgumentParser(description="Summarize peak counts and widths.")
    p.add_argument("--peaks", nargs="+", required=True, help="Gzipped peak files.")
    p.add_argument("--labels", nargs="+", required=True, help="One label per file.")
    p.add_argument("--title", default="", help="Figure title.")
    p.add_argument("--table", required=True, help="TSV to write.")
    p.add_argument("--plot", required=True, help="SVG to write (PNG beside it).")
    return p


def widths(path):
    """Read the region widths of a peak file.

    Parameters
    ----------
    path : str
        Gzipped peak file, possibly empty.

    Returns
    -------
    polars.Series
        End minus start of every region.
    """
    try:
        peaks = pl.read_csv(path, separator="\t", has_header=False, columns=[1, 2])
    except pl.exceptions.NoDataError:
        return pl.Series("width", [], dtype=pl.Int64)
    return (peaks[:, 1] - peaks[:, 0]).alias("width")


def summarize(label, width):
    """Summarize one peak set.

    Parameters
    ----------
    label : str
        Peak set label.
    width : polars.Series
        Region widths.

    Returns
    -------
    dict of str to object
        Label, peak count and the `STATS` columns, null for an empty set.
    """
    return {
        "label": label,
        "num_peaks": width.len(),
        "min_size": width.min(),
        "25_pct": width.quantile(0.25, "linear"),
        "50_pct": width.quantile(0.5, "linear"),
        "75_pct": width.quantile(0.75, "linear"),
        "max_size": width.max(),
        "mean": width.mean(),
    }


def main(argv=None):
    """Write the table and draw each set's width distribution.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    sets = {label: widths(p) for label, p in zip(args.labels, args.peaks, strict=True)}
    pl.DataFrame(
        [summarize(label, w) for label, w in sets.items()],
        schema_overrides={s: pl.Float64 for s in STATS},
    ).write_csv(args.table, separator="\t")

    pk.use("plain")
    fig, ax = plt.subplots(figsize=(6, 4))
    for label, width in sets.items():
        if width.len():
            counts, edges = np.histogram(width.to_numpy(), bins=100)
            ax.plot((edges[1:] + edges[:-1]) / 2, counts, label=label)
    ax.set(xlabel="Peak width (bp)", ylabel="Peaks", title=args.title)
    ax.legend()
    pk.despine(ax)
    pk.save_figure(fig, args.plot)
    plt.close(fig)


if __name__ == "__main__":
    main()
