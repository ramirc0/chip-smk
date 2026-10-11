"""Gather a run's QC into one long table and one figure.

Rows are ``experiment, section, unit, metric, value``. ``section.unit.metric``
follows the paths of ENCODE's ``qc.json``, so the table lines up with ENCODE's
published QC. File names are ``<exp>.<unit>...``; names never contain dots.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import orjson
import poikilos as pk
import polars as pl

# samtools flagstat JSON fields, as ENCODE's samstat keys.
FLAGSTAT = {
    "total": "total_reads",
    "mapped": "mapped_reads",
    "mapped %": "pct_mapped_reads",
    "paired in sequencing": "paired_reads",
    "read1": "read1",
    "read2": "read2",
    "properly paired": "properly_paired_reads",
    "properly paired %": "pct_properly_paired_reads",
    "with itself and mate mapped": "with_itself",
    "singletons": "singletons",
    "singletons %": "pct_singletons",
    "with mate mapped to a different chr": "diff_chroms",
}
XCOR = [
    "subsampled_reads",
    "estimated_fragment_len",
    "corr_estimated_fragment_len",
    "phantom_peak",
    "corr_phantom_peak",
    "argmin_corr",
    "min_corr",
    "NSC",
    "RSC",
]
JSD = [
    "auc",
    "syn_auc",
    "x_intercept",
    "syn_x_intercept",
    "elbow_pt",
    "syn_elbow_pt",
    "jsd",
    "syn_jsd",
    "pct_genome_enrich",
    "diff_enrich",
    "ch_div",
]
# Panels of the summary figure: section, metric, title.
PANELS = [
    ("align.samstat", "pct_mapped_reads", "Mapped reads (%)"),
    ("align.dup", "pct_duplicate_reads", "Duplicate reads (%)"),
    ("lib_complexity.lib_complexity", "NRF", "NRF"),
    ("lib_complexity.lib_complexity", "PBC1", "PBC1"),
    ("align_enrich.xcor_score", "NSC", "NSC"),
    ("align_enrich.xcor_score", "RSC", "RSC"),
    ("align_enrich.jsd", "syn_jsd", "Synthetic JSD"),
    ("replication.reproducibility", "N_opt", "Optimal peaks"),
]


def build_parser():
    """Build the qc_summary.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the qc_summary.py flags.
    """
    p = argparse.ArgumentParser(description="Gather a run's QC.")
    for flag, what in (
        ("--flagstat", "samtools flagstat JSON (align/, filter/, .dup.)"),
        ("--pbc", "pbc_qc.py JSON"),
        ("--fastp", "fastp JSON, one per lane"),
        ("--xcor", "run_spp.R cc.qc"),
        ("--jsd", "plotFingerprint quality metrics"),
        ("--frip", "frip TSV"),
        ("--reproducibility", "peak_sets.py TSV"),
        ("--peak-qc", "peak_qc.py TSV"),
    ):
        p.add_argument(flag, nargs="*", default=[], help=f"{what} files.")
    p.add_argument("--table", required=True, help="TSV to write.")
    p.add_argument("--plot", required=True, help="SVG to write (PNG beside it).")
    return p


def name_parts(path):
    """Split a QC file name into experiment and unit.

    Parameters
    ----------
    path : str
        File named ``<exp>.<unit>.<rest>``.

    Returns
    -------
    tuple of str
        Experiment and unit.
    """
    exp, unit = Path(path).name.split(".")[:2]
    return exp, unit


def flagstat_rows(path):
    """Read a flagstat JSON as samstat or duplicate rows.

    Parameters
    ----------
    path : str
        ``align/<exp>.<lib>.flagstat.json``, ``filter/<exp>.<lib>.flagstat.json``
        or ``filter/<exp>.<lib>.dup.flagstat.json``.

    Returns
    -------
    list of tuple
        ``(experiment, section, unit, metric, value)`` rows.
    """
    exp, lib = name_parts(path)
    with open(path, "rb") as f:
        stats = orjson.loads(f.read())["QC-passed reads"]
    if path.endswith(".dup.flagstat.json"):
        # Picard's MarkDuplicates fields, from the duplicate-marked reads.
        total, dups, paired = (
            stats["total"],
            stats["duplicates"],
            stats["paired in sequencing"],
        )
        metrics = {
            "unpaired_reads": total - paired,
            "paired_reads": paired // 2,
            "unpaired_duplicate_reads": dups if not paired else 0,
            "paired_duplicate_reads": dups // 2 if paired else 0,
            "pct_duplicate_reads": 100 * dups / total if total else 0.0,
        }
        section = "align.dup"
    else:
        metrics = {key: stats[field] for field, key in FLAGSTAT.items()}
        section = (
            "align.samstat"
            if Path(path).parent.name == "align"
            else "align.nodup_samstat"
        )
    return [(exp, section, lib, k, v) for k, v in metrics.items()]


def pbc_rows(path):
    """Read a library complexity JSON.

    Parameters
    ----------
    path : str
        ``pbc_qc/<exp>.<lib>.json``.

    Returns
    -------
    list of tuple
        ``(experiment, section, unit, metric, value)`` rows.
    """
    exp, lib = name_parts(path)
    with open(path, "rb") as f:
        metrics = orjson.loads(f.read())
    return [
        (exp, "lib_complexity.lib_complexity", lib, k, v) for k, v in metrics.items()
    ]


def fastp_rows(paths):
    """Sum the fastp reports of each library's lanes.

    Parameters
    ----------
    paths : list of str
        ``fastp/<exp>.<lib>.<lane>.json``.

    Returns
    -------
    list of tuple
        ``(experiment, section, unit, metric, value)`` rows: total reads and
        Q30 base fraction per library.
    """
    totals = {}
    for path in paths:
        with open(path, "rb") as f:
            before = orjson.loads(f.read())["summary"]["before_filtering"]
        t = totals.setdefault(name_parts(path), [0, 0, 0])
        t[0] += before["total_reads"]
        t[1] += before["total_bases"]
        t[2] += before["q30_bases"]
    rows = []
    for (exp, lib), (reads, bases, q30) in totals.items():
        rows.append((exp, "fastp", lib, "total_reads", reads))
        rows.append((exp, "fastp", lib, "q30_rate", q30 / bases if bases else 0.0))
    return rows


def xcor_rows(path):
    """Read a run_spp.R cross-correlation row.

    Parameters
    ----------
    path : str
        ``xcor/<exp>.<rep>.cc.qc``.

    Returns
    -------
    list of tuple
        ``(experiment, section, unit, metric, value)`` rows.
    """
    exp, rep = name_parts(path)
    with open(path) as f:
        fields = f.read().rstrip("\n").split("\t")[1:]
    return [(exp, "align_enrich.xcor_score", rep, k, v) for k, v in zip(XCOR, fields)]


def jsd_rows(path):
    """Read plotFingerprint quality metrics of the ChIP replicates.

    Parameters
    ----------
    path : str
        ``jsd/<exp>.tsv``. Control rows are left out, as ENCODE does.

    Returns
    -------
    list of tuple
        ``(experiment, section, unit, metric, value)`` rows.
    """
    exp = Path(path).name.split(".")[0]
    rows = []
    with open(path) as f:
        for line in f.read().splitlines()[1:]:
            label, *values = line.split("\t")
            if label.startswith("rep"):
                rows += [
                    (exp, "align_enrich.jsd", label, k, v) for k, v in zip(JSD, values)
                ]
    return rows


def frip_rows(path):
    """Read a FRiP table as rows under the caller, overlap or IDR.

    Parameters
    ----------
    path : str
        ``frip/<exp>.<sample>.tsv`` with ``peak_set`` labels ``<caller>`` or
        ``<method>.<pair>``.

    Returns
    -------
    list of tuple
        ``(experiment, section, unit, metric, value)`` rows.
    """
    exp, sample = name_parts(path)
    rows = []
    for r in pl.read_csv(path, separator="\t").iter_rows(named=True):
        method, _, pair = r["peak_set"].partition(".")
        frip = r["reads_in_peaks"] / r["reads"] if r["reads"] else 0.0
        section = f"peak_enrich.frac_reads_in_peaks.{method}"
        rows.append((exp, section, pair or sample, "frip", frip))
    return rows


def reproducibility_rows(path):
    """Read a reproducibility table.

    Parameters
    ----------
    path : str
        ``reproducibility/<exp>.<method>.<type>.tsv``.

    Returns
    -------
    list of tuple
        ``(experiment, section, unit, metric, value)`` rows.
    """
    exp, method = name_parts(path)
    row = pl.read_csv(path, separator="\t").row(0, named=True)
    return [(exp, "replication.reproducibility", method, k, v) for k, v in row.items()]


def peak_qc_rows(path):
    """Read a peak QC table as peak counts and region size statistics.

    Parameters
    ----------
    path : str
        ``peak_qc/<exp>.tsv``.

    Returns
    -------
    list of tuple
        ``(experiment, section, unit, metric, value)`` rows.
    """
    exp = Path(path).name.split(".")[0]
    rows = []
    for r in pl.read_csv(path, separator="\t").iter_rows(named=True):
        label = r.pop("label")
        rows.append(
            (exp, "replication.num_peaks", label, "num_peaks", r.pop("num_peaks"))
        )
        rows += [(exp, "peak_stat.peak_region_size", label, k, v) for k, v in r.items()]
    return rows


def main(argv=None):
    """Gather every QC file into the table, then draw the summary panels.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    rows = fastp_rows(args.fastp)
    for paths, read in (
        (args.flagstat, flagstat_rows),
        (args.pbc, pbc_rows),
        (args.xcor, xcor_rows),
        (args.jsd, jsd_rows),
        (args.frip, frip_rows),
        (args.reproducibility, reproducibility_rows),
        (args.peak_qc, peak_qc_rows),
    ):
        for path in paths:
            rows += read(path)
    table = pl.DataFrame(
        [(e, s, u, m, "" if v is None else str(v)) for e, s, u, m, v in rows],
        schema=["experiment", "section", "unit", "metric", "value"],
        orient="row",
    )
    table.write_csv(args.table, separator="\t")
    _plot(table, args.plot)


def _plot(table, path):
    """Draw one bar panel per summary metric, a bar per experiment and unit.

    Parameters
    ----------
    table : polars.DataFrame
        The QC table.
    path : str
        SVG to write. The PNG goes beside it.
    """
    pk.use("plain")
    fig, axes = plt.subplots(2, len(PANELS) // 2, figsize=(4 * len(PANELS) // 2, 8))
    for ax, (section, metric, title) in zip(axes.flat, PANELS):
        values = table.filter(
            (pl.col("section") == section) & (pl.col("metric") == metric)
        ).with_columns(pl.col("value").cast(pl.Float64, strict=False))
        labels = [f"{e} {u}" for e, u in zip(values["experiment"], values["unit"])]
        ax.barh(labels, values["value"].fill_null(0).to_list())
        ax.set(title=title)
        ax.invert_yaxis()
        pk.despine(ax, categorical_y=True)
    pk.save_figure(fig, path)
    plt.close(fig)


if __name__ == "__main__":
    main()
