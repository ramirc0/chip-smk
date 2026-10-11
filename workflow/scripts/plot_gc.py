"""Plot Picard's GC bias metrics: normalized coverage and base quality by GC."""

import argparse
import io

import matplotlib.pyplot as plt
import poikilos as pk
import polars as pl


def build_parser():
    """Build the plot_gc.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_gc.py flags.
    """
    p = argparse.ArgumentParser(description="Plot GC bias.")
    p.add_argument(
        "--metrics", required=True, help="Picard CollectGcBiasMetrics output."
    )
    p.add_argument("--plot", required=True, help="SVG to write (PNG beside it).")
    return p


def main(argv=None):
    """Draw ENCODE's plot_gc series (coverage, genome windows, base quality) by GC.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    # Picard puts a blank line between its header comments and the table.
    with open(args.metrics) as f:
        table = "".join(line for line in f if line.strip() and line[0] != "#")
    gc = pl.read_csv(io.StringIO(table), separator="\t")
    pk.use("plain")
    # Stacked panels instead of ENCODE's twin y axes: the plain look has no right spine.
    fig, axes = plt.subplots(3, 1, figsize=(7, 8), sharex=True, height_ratios=(2, 1, 1))
    for ax, column, label in zip(
        axes,
        ("NORMALIZED_COVERAGE", "WINDOWS", "MEAN_BASE_QUALITY"),
        ("Normalized coverage", "Genome windows (fraction)", "Mean base quality"),
    ):
        y = gc[column] / gc[column].sum() if column == "WINDOWS" else gc[column]
        ax.plot(gc["GC"], y)
        ax.set_ylabel(label)
        pk.despine(ax)
    axes[0].set_title("GC bias")
    axes[-1].set_xlabel("GC (%)")
    pk.save_figure(fig, args.plot)
    plt.close(fig)


if __name__ == "__main__":
    main()
