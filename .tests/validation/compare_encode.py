"""Print the validation run's QC beside ENCODE's published values.

ENCODE ran older tools (bowtie2 2.3.4.3, MACS2 2.2.4, SPP 1.15) and marked
duplicates with Picard, so values are close, not equal. The integration tests
compare step by step on identical inputs. Run from the repository root after the
validation run:

    pixi run --manifest-path workflow/envs/pixi.toml python .tests/validation/compare_encode.py
"""

import argparse
import math
import subprocess
from pathlib import Path

import polars as pl
import pyBigWig
from zlib_ng import gzip_ng

ENCODE_DATA = Path("/data/projects/encode/data")


def build_parser():
    """Build the compare_encode.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the compare_encode.py flags.
    """
    p = argparse.ArgumentParser(description="Compare QC with ENCODE's.")
    p.add_argument("--run", default="results/encode-validation", help="Run directory.")
    p.add_argument(
        "--encode", default=".tests/validation/encode_qc.tsv", help="ENCODE QC."
    )
    return p


def peak_count(path):
    """Count the peaks in a gzipped peak file.

    Parameters
    ----------
    path : pathlib.Path
        Gzipped peak file.

    Returns
    -------
    int
        Lines in the file.
    """
    with gzip_ng.open(path, "rb") as f:
        return sum(1 for _ in f)


def peaks_hit(a, b):
    """Return the fraction of peaks in ``a`` that overlap a peak in ``b``.

    Parameters
    ----------
    a, b : pathlib.Path
        Gzipped peak files.

    Returns
    -------
    float
        Overlapping peaks over all peaks of ``a``.
    """
    hits = subprocess.run(
        ["bedtools", "intersect", "-u", "-a", str(a), "-b", str(b)],
        capture_output=True,
        text=True,
        check=True,
    )
    return len(hits.stdout.splitlines()) / peak_count(a)


def bigwig_stats(path):
    """Return a bigWig's mean, max and standard deviation, as ``bigWigInfo`` reports them.

    Parameters
    ----------
    path : pathlib.Path
        bigWig file.

    Returns
    -------
    dict of str to float
        ``mean``, ``max`` and ``std`` over the covered bases, from the header's sums.
    """
    bw = pyBigWig.open(str(path))
    h = bw.header()
    bw.close()
    n, total = h["nBasesCovered"], h["sumData"]
    return {
        "mean": total / n,
        "max": h["maxVal"],
        "std": math.sqrt((h["sumSquared"] - total * total / n) / (n - 1)),
    }


def our_value(ours, row, sets):
    """Look up our value for one ENCODE row.

    Parameters
    ----------
    ours : polars.DataFrame
        The run's ``qc.tsv``.
    row : dict of str to str
        ENCODE row: ``experiment, section, unit, metric, value``.
    sets : dict of str to str
        Our pair label of the ``optimal`` and ``conservative`` sets, per
        experiment and method.

    Returns
    -------
    str
        Our value, or ``-`` when the run has none.
    """
    unit = row["unit"]
    if row["section"].startswith("peak_enrich."):
        method = row["section"].rsplit(".", 1)[1]
        unit = sets.get((row["experiment"], method, unit), unit)
    match = ours.filter(
        (pl.col("experiment") == row["experiment"])
        & (pl.col("section") == row["section"])
        & (pl.col("unit") == unit)
        & (pl.col("metric") == row["metric"])
    )["value"]
    return match[0] if match.len() else "-"


def main(argv=None):
    """Print each metric per experiment, ours beside ENCODE's, then peak agreement.

    The pooled signal tracks follow, as their mean, max and standard deviation.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    run = Path(args.run)
    ours = pl.read_csv(run / "qc.tsv", separator="\t", infer_schema=False)
    encode = pl.read_csv(
        args.encode, separator="\t", comment_prefix="#", infer_schema=False
    )
    chosen = ours.filter(
        (pl.col("section") == "replication.reproducibility")
        & pl.col("metric").is_in(["opt_set", "consv_set"])
    )
    sets = {
        (e, u, "optimal" if m == "opt_set" else "conservative"): v
        for e, u, m, v in chosen.select("experiment", "unit", "metric", "value").rows()
    }
    methods = dict(
        encode.filter(pl.col("section") == "replication.reproducibility")
        .select("experiment", "unit")
        .unique()
        .rows()
    )

    for (exp,), rows in encode.group_by("experiment", maintain_order=True):
        print(f"\n{exp}\n{'metric':<48}{'ours':>14}{'ENCODE':>14}")
        metrics = rows.filter(~pl.col("section").is_in(["files", "signal"]))
        for row in metrics.iter_rows(named=True):
            label = f"{row['section'].split('.')[-1]} {row['unit']} {row['metric']}"
            mine = our_value(ours, row, sets)
            print(f"{label:<48}{_short(mine):>14}{_short(row['value']):>14}")
        method = methods[exp]
        for kind, path in (
            rows.filter(pl.col("section") == "files").select("unit", "value").rows()
        ):
            mine = next(
                (run / exp / "reproducibility").glob(f"{exp}.{method}.{kind}.*.gz")
            )
            theirs = ENCODE_DATA / path
            print(
                f"{kind} peaks: ours in ENCODE's {peaks_hit(mine, theirs):.3f}, "
                f"ENCODE's in ours {peaks_hit(theirs, mine):.3f}"
            )
        signal = rows.filter(pl.col("section") == "signal")
        for unit, track, path in signal.select("unit", "metric", "value").rows():
            mine = bigwig_stats(run / exp / "signal" / f"{exp}.{unit}.{track}.bigWig")
            theirs = bigwig_stats(ENCODE_DATA / path)
            for stat in ("mean", "max", "std"):
                label = f"signal {unit} {track} {stat}"
                print(f"{label:<48}{mine[stat]:>14.4g}{theirs[stat]:>14.4g}")


def _short(value):
    """Format a value for the table: numbers to four significant digits.

    Parameters
    ----------
    value : str
        Value as text.

    Returns
    -------
    str
        The formatted value.
    """
    try:
        return f"{float(value):.4g}"
    except ValueError:
        return value


if __name__ == "__main__":
    main()
