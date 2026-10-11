"""Aggregate per-sample metrics into a distribution plot, one facet per metric."""

import argparse

import matplotlib.pyplot as plt
import poikilos as pk
import polars as pl


def build_parser():
    """Build the plot_summary.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_summary.py flags.
    """
    p = argparse.ArgumentParser(description="Plot the metric distribution.")
    p.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Per-sample metrics TSVs.",
    )
    p.add_argument(
        "-o", "--output", required=True, help="Summary SVG to write (PNG beside it)."
    )
    return p


def main(argv=None):
    """Concatenate the metrics tables and plot one distribution per metric.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    table = pl.concat([pl.read_csv(f, separator="\t") for f in args.input])
    metrics = table["metric"].unique().sort().to_list()

    pk.use("plain")
    fig, axes = plt.subplots(
        1, len(metrics), figsize=(3 * len(metrics), 3), squeeze=False
    )
    for ax, metric in zip(axes[0], metrics):
        vals = table.filter(pl.col("metric") == metric)["value"].to_list()
        ax.boxplot(vals)
        ax.set(title=metric, ylabel="value", xticks=[])
        pk.despine(ax, categorical_x=True)
    fig.suptitle("Metric distribution across samples")
    pk.save_figure(fig, args.output)
    plt.close(fig)


if __name__ == "__main__":
    main()
