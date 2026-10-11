"""Plot a strand cross-correlation profile from run_spp.R.

Marks the estimated fragment length and the phantom peak (read length) from the
``cc.qc`` row, as run_spp.R's own plot does.
"""

import argparse

import matplotlib.pyplot as plt
import poikilos as pk
import polars as pl


def build_parser():
    """Build the plot_xcor.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_xcor.py flags.
    """
    p = argparse.ArgumentParser(description="Plot a cross-correlation profile.")
    p.add_argument("--qc", required=True, help="run_spp.R -out row, first estimates.")
    p.add_argument(
        "--profile", required=True, help="Profile TSV: shift x, correlation y."
    )
    p.add_argument("--title", default="", help="Figure title.")
    p.add_argument("--output", required=True, help="SVG to write (PNG beside it).")
    return p


def main(argv=None):
    """Draw the profile with the fragment length and phantom peak marked.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    profile = pl.read_csv(args.profile, separator="\t")
    with open(args.qc) as f:
        qc = f.read().split("\t")
    fraglen, phantom, nsc, rsc = int(qc[2]), int(qc[4]), float(qc[8]), float(qc[9])

    pk.use("plain")
    colors = pk.PALETTES["plain"].colors
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(profile["x"], profile["y"])
    ax.axvline(
        fraglen, color=colors["red"], linestyle="--", label=f"fragment {fraglen}"
    )
    ax.axvline(phantom, color=colors["blue"], linestyle=":", label=f"phantom {phantom}")
    ax.set(
        xlabel="Strand shift (bp)",
        ylabel="Cross-correlation",
        title=f"{args.title}\nNSC {nsc:.3f}, RSC {rsc:.3f}",
    )
    ax.legend()
    pk.despine(ax)
    pk.save_figure(fig, args.output)
    plt.close(fig)


if __name__ == "__main__":
    main()
