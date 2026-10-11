"""Pick an experiment's conservative and optimal peak sets, as ENCODE's reproducibility does.

Takes the filtered overlap or IDR peaks of each comparison, named
``<exp>.<pair>.bfilt.<type>.gz``. Writes ENCODE's ``reproducibility.qc``
columns, with infinite ratios instead of errors for empty peak sets.
"""

import argparse
import math
import shutil
from pathlib import Path

import polars as pl
from zlib_ng import gzip_ng


def build_parser():
    """Build the peak_sets.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the peak_sets.py flags.
    """
    p = argparse.ArgumentParser(description="Pick conservative and optimal peak sets.")
    p.add_argument(
        "--true-pairs",
        nargs="*",
        default=[],
        help="Peaks of each true replicate pair, in ENCODE's order. Empty for one replicate.",
    )
    p.add_argument(
        "--self-pseudo",
        nargs="+",
        required=True,
        help="Peaks of each replicate's pseudoreplicate pair.",
    )
    p.add_argument(
        "--pooled-pseudo",
        nargs="*",
        default=[],
        help="Peaks of the pooled pseudoreplicate pair. Empty for one replicate.",
    )
    p.add_argument("--conservative", required=True, help="Conservative set to write.")
    p.add_argument("--optimal", required=True, help="Optimal set to write.")
    p.add_argument("--table", required=True, help="Reproducibility TSV to write.")
    return p


def main(argv=None):
    """Count peaks, pick the two sets, and write them with the QC table.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    n_self = [_count(p) for p in args.self_pseudo]
    if args.true_pairs:
        n_true = [_count(p) for p in args.true_pairs]
        nt, np_ = max(n_true), _count(args.pooled_pseudo[0])
        conservative = args.true_pairs[n_true.index(nt)]
        optimal = conservative if nt > np_ else args.pooled_pseudo[0]
        rescue = _ratio(nt, np_)
        self_consistency = _ratio(max(n_self), min(n_self))
    else:
        nt = np_ = 0
        conservative = optimal = args.self_pseudo[0]
        rescue, self_consistency = 0.0, 1.0

    shutil.copyfile(conservative, args.conservative)
    shutil.copyfile(optimal, args.optimal)

    failed = (rescue > 2) + (self_consistency > 2)
    row = {"Nt": nt}
    row |= {f"N{i}": n for i, n in enumerate(n_self, start=1)}
    row |= {
        "Np": np_,
        "N_opt": _count(optimal),
        "N_consv": _count(conservative),
        "opt_set": _label(optimal),
        "consv_set": _label(conservative),
        "rescue_ratio": rescue,
        "self_consistency_ratio": self_consistency,
        "reproducibility": ["pass", "borderline", "fail"][failed],
    }
    pl.DataFrame([row]).write_csv(args.table, separator="\t")


def _count(path):
    """Count the peaks in a gzipped peak file.

    Parameters
    ----------
    path : str
        Gzipped peak file.

    Returns
    -------
    int
        Number of lines.
    """
    with gzip_ng.open(path, "rb") as f:
        return sum(1 for _ in f)


def _label(path):
    """Return the pair label of a peak file.

    Parameters
    ----------
    path : str
        File named ``<exp>.<pair>.bfilt.<type>.gz``.

    Returns
    -------
    str
        The pair, such as ``rep1_vs_rep2`` or ``pooled-pr1_vs_pooled-pr2``.
    """
    return Path(path).name.split(".")[1]


def _ratio(a, b):
    """Return the larger count over the smaller.

    Parameters
    ----------
    a, b : int
        Peak counts.

    Returns
    -------
    float
        `max / min`, or infinity when the smaller count is zero.
    """
    lo, hi = sorted((a, b))
    return hi / lo if lo else math.inf


if __name__ == "__main__":
    main()
