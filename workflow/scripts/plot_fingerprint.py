"""Plot fingerprints from deeptools plotFingerprint's raw bin counts.

Each sample's bins are sorted by read count. Its curve is the fraction of reads
in the lowest-count fraction of bins. Uniform coverage follows the diagonal.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import poikilos as pk
import polars as pl


def build_parser():
    """Build the plot_fingerprint.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_fingerprint.py flags.
    """
    p = argparse.ArgumentParser(description="Plot fingerprints.")
    p.add_argument("--counts", required=True, help="plotFingerprint --outRawCounts.")
    p.add_argument("--title", default="", help="Figure title.")
    p.add_argument("--plot", required=True, help="SVG to write (PNG beside it).")
    return p


def main(argv=None):
    """Draw one cumulative read fraction curve per sample.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    raw = pl.read_csv(args.counts, separator="\t", comment_prefix="#", quote_char="'")
    pk.use("plain")
    fig, ax = plt.subplots(figsize=(5, 5))
    for label in raw.columns:
        counts = np.sort(raw[label].to_numpy())
        ax.plot(
            np.linspace(0, 1, counts.size),
            np.cumsum(counts) / counts.sum(),
            label=label,
        )
    ax.plot([0, 1], [0, 1], color=pk.PALETTES["plain"].foreground, linestyle=":")
    ax.set(xlabel="Fraction of bins", ylabel="Fraction of reads", title=args.title)
    ax.legend()
    pk.despine(ax)
    pk.save_figure(fig, args.plot)
    plt.close(fig)


if __name__ == "__main__":
    main()
