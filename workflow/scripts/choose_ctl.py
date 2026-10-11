"""Choose each replicate's control and its subsample depth, as ENCODE's choose_ctl does.

Reads depths from tagAlign count files (``lines<TAB>bytes``) and writes one row per
replicate and one for the pool: the control to call against (``ctlN`` or
``pooled``) and the reads to subsample it to (0 = all).
"""

import argparse

import polars as pl


def build_parser():
    """Build the choose_ctl.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the choose_ctl.py flags.
    """
    p = argparse.ArgumentParser(description="Choose controls for ChIP replicates.")
    p.add_argument(
        "--reps", nargs="+", required=True, help="Replicate tagAlign count files."
    )
    p.add_argument(
        "--controls", nargs="+", required=True, help="Control tagAlign count files."
    )
    p.add_argument(
        "--ctl-depth-ratio",
        type=float,
        default=1.2,
        help="Controls differing in depth by more than this are pooled.",
    )
    p.add_argument(
        "--ctl-depth-limit",
        type=int,
        default=200_000_000,
        help="Control reads kept at most, unless the ratio limit allows more.",
    )
    p.add_argument(
        "--exp-ctl-depth-ratio-limit",
        type=float,
        default=5.0,
        help="Control reads kept at most, as a multiple of the replicate's.",
    )
    p.add_argument(
        "--always-use-pooled-ctl",
        action="store_true",
        help="Call every replicate against the pooled control.",
    )
    p.add_argument(
        "--paired-end-control",
        action="store_true",
        help="Round subsample depths down to whole pairs.",
    )
    p.add_argument("--output", required=True, help="TSV to write.")
    return p


def choose(depths, ctl_depths, ratio, always_pooled):
    """Pick a control for each replicate.

    Parameters
    ----------
    depths : list of int
        Reads per replicate.
    ctl_depths : list of int
        Reads per control.
    ratio : float
        Maximum depth ratio between two controls before all are pooled.
    always_pooled : bool
        Call every replicate against the pooled control.

    Returns
    -------
    list of int
        Control index per replicate; -1 is the pooled control.
    """
    n = len(depths)
    if len(ctl_depths) == 1:
        return [0] * n
    if always_pooled:
        return [-1] * n
    uneven = any(
        a / b > ratio or b / a > ratio
        for i, a in enumerate(ctl_depths)
        for b in ctl_depths[i + 1 :]
    )
    if uneven:
        return [-1] * n
    return [
        i if i < len(ctl_depths) and ctl_depths[i] >= depths[i] else -1
        for i in range(n)
    ]


def subsample(depth, ctl_depth, ratio_limit, depth_limit, paired):
    """Return the reads to subsample a control to, or 0 to keep it whole.

    Parameters
    ----------
    depth : int
        Reads in the ChIP tagAlign.
    ctl_depth : int
        Reads in the chosen control.
    ratio_limit : float
        Control reads allowed per ChIP read.
    depth_limit : int
        Control reads always allowed.
    paired : bool
        Round down to an even number of reads (whole pairs).

    Returns
    -------
    int
        Subsample size, or 0.
    """
    if not (ratio_limit or depth_limit):
        return 0
    limit = int(max(depth * ratio_limit, depth_limit))
    if ctl_depth <= limit:
        return 0
    return limit - limit % 2 if paired else limit


def main(argv=None):
    """Read the depths, choose controls and write the table.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    depths = [_lines(p) for p in args.reps]
    ctl_depths = [_lines(p) for p in args.controls]
    chosen = choose(
        depths, ctl_depths, args.ctl_depth_ratio, args.always_use_pooled_ctl
    )

    rows = []
    for i, c in enumerate(chosen):
        ctl_depth = ctl_depths[c] if c >= 0 else sum(ctl_depths)
        rows.append(
            (f"rep{i + 1}", "pooled" if c < 0 else f"ctl{c + 1}", ctl_depth, depths[i])
        )
    # The pooled replicate uses the only control, else the pooled one.
    pooled_ctl = "ctl1" if len(ctl_depths) == 1 else "pooled"
    rows.append(("pooled", pooled_ctl, sum(ctl_depths), sum(depths)))

    pl.DataFrame(
        {
            "sample": [r[0] for r in rows],
            "control": [r[1] for r in rows],
            "subsample": [
                subsample(
                    depth,
                    ctl_depth,
                    args.exp_ctl_depth_ratio_limit,
                    args.ctl_depth_limit,
                    args.paired_end_control,
                )
                for _, _, ctl_depth, depth in rows
            ],
        }
    ).write_csv(args.output, separator="\t")


def _lines(path):
    """Read the line count of a tagAlign count file.

    Parameters
    ----------
    path : str
        File holding ``lines<TAB>bytes``.

    Returns
    -------
    int
        Lines in the tagAlign.
    """
    with open(path) as f:
        return int(f.read().split()[0])


if __name__ == "__main__":
    main()
