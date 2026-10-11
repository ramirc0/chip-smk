"""Per-sample analysis: read a value column, write metrics + a curve plot.

Example script. Keep the shape: a `build_parser()` returning the parser and a
`main()` reading the declared inputs / writing the declared outputs.
"""

import argparse

import matplotlib.pyplot as plt
import poikilos as pk
import polars as pl


def build_parser():
    """Build the analyze.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the analyze.py flags.
    """
    p = argparse.ArgumentParser(description="Compute per-sample metrics.")
    p.add_argument("-i", "--input", required=True, help="Prepared input TSV.")
    p.add_argument("-o", "--output", required=True, help="Metrics TSV to write.")
    p.add_argument(
        "--plot", required=True, help="Per-sample curve SVG to write (PNG beside it)."
    )
    p.add_argument("--sample", required=True, help="Sample id (used as plot title).")
    p.add_argument("--window", type=int, default=100, help="Smoothing window.")
    p.add_argument("--threshold", type=float, default=0.5, help="Call threshold.")
    p.add_argument(
        "--metric",
        default="mean",
        choices=["mean", "median"],
        help="Summary statistic.",
    )
    p.add_argument("--normalize", action="store_true", help="Min-max normalize.")
    p.add_argument("--random_state", type=int, default=None, help="RNG seed.")
    return p


def main(argv=None):
    """Read the input value column, compute metrics, write the table + plot.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    values = pl.read_csv(args.input, separator="\t")["value"]
    if args.normalize:
        lo, hi = values.min(), values.max()
        values = (values - lo) / (hi - lo) if hi > lo else values * 0

    stat = values.mean() if args.metric == "mean" else values.median()
    pl.DataFrame(
        {
            "sample": [args.sample],
            "metric": [args.metric],
            "value": [stat],
            "n_above_threshold": [int((values > args.threshold).sum())],
            "count": [values.len()],
        }
    ).write_csv(args.output, separator="\t")

    _plot_curve(values.sort(), args.sample, args.plot)


def _plot_curve(values, sample, path):
    """Plot sorted values against their rank.

    Parameters
    ----------
    values : polars.Series
        Values in ascending order.
    sample : str
        Sample ID, used as the title.
    path : str
        Curve SVG to write. The PNG goes beside it.
    """
    pk.use("plain")
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot(range(values.len()), values.to_list())
    ax.set(xlabel="rank", ylabel="value", title=sample)
    pk.despine(ax)
    pk.save_figure(fig, path)
    plt.close(fig)


if __name__ == "__main__":
    main()
