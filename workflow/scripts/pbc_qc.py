"""Compute library complexity (NRF, PBC1, PBC2) from bedtools bamtobed, as ENCODE does."""

import argparse
import sys

import orjson
import polars as pl

# bedtools bamtobed columns, and the ENCODE keys a read or fragment is counted on.
BED = ["chrom", "start", "end", "name", "score", "strand"]
BEDPE = [
    "chrom1",
    "start1",
    "end1",
    "chrom2",
    "start2",
    "end2",
    "name",
    "score",
    "strand1",
    "strand2",
]
KEY = {
    "se": ["chrom", "start", "end", "strand"],
    "pe": ["chrom1", "start1", "chrom2", "end2", "strand1", "strand2"],
}


def build_parser():
    """Build the pbc_qc.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the pbc_qc.py flags.
    """
    p = argparse.ArgumentParser(description="Compute NRF, PBC1 and PBC2.")
    p.add_argument(
        "--paired-end",
        action="store_true",
        help="Read `bamtobed -bedpe` fragments instead of single reads.",
    )
    p.add_argument("--mito", default="chrM", help="Mitochondrial contig to skip.")
    p.add_argument("--output", required=True, help="JSON to write.")
    return p


def library_complexity(reads, key, mito):
    """Count reads by position and derive ENCODE's complexity metrics.

    Parameters
    ----------
    reads : polars.DataFrame
        One row per read or fragment, with the `key` columns.
    key : list of str
        Columns a position is keyed on. The first is the contig; rows on
        ``mito`` are skipped.
    mito : str
        Mitochondrial contig.

    Returns
    -------
    dict of str to float
        Totals and NRF, PBC1 and PBC2. PBC2 is -1 when no position is seen
        exactly twice, and a ratio is 0 when its denominator is 0.
    """
    counts = reads.filter(pl.col(key[0]) != mito).group_by(key).len()["len"]
    total = int(counts.sum())
    distinct = counts.len()
    one = int((counts == 1).sum())
    two = int((counts == 2).sum())
    return {
        "total_fragments": total,
        "distinct_fragments": distinct,
        "positions_with_one_read": one,
        "positions_with_two_reads": two,
        "NRF": distinct / total if total else 0.0,
        "PBC1": one / distinct if distinct else 0.0,
        "PBC2": one / two if two else -1.0,
    }


def main(argv=None):
    """Read bamtobed output from stdin and write the complexity metrics as JSON.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    layout = "pe" if args.paired_end else "se"
    columns = BEDPE if args.paired_end else BED
    reads = pl.read_csv(
        sys.stdin.buffer,
        separator="\t",
        has_header=False,
        new_columns=columns,
        schema_overrides={c: pl.String for c in columns if c.startswith("chrom")},
    ).select(KEY[layout])
    metrics = library_complexity(reads, KEY[layout], args.mito)
    with open(args.output, "wb") as f:
        f.write(orjson.dumps(metrics, option=orjson.OPT_INDENT_2))


if __name__ == "__main__":
    main()
