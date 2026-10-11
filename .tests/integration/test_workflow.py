"""Check sheet validation, failure behavior and the shape of the job graph."""

import gzip
import re

import polars as pl
import pytest
from conftest import CHRSZ, CONFIG, ROOT, lines, run_rule, snakemake

HEADER = "experiment\treplicate\ttype\tcontrol\tcrop_length\tr1\tr2\n"
CHIP = "Chip\t1\ttf\tInput\t\tchip_R1.fastq.gz\t\n"
INPUT = "Input\t1\tcontrol\t\t\tinput_R1.fastq.gz\t\n"


def dry_run(tmp_path, sheet, *args, config=""):
    """Dry-run the workflow on a sheet in a scratch directory.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Scratch directory.
    sheet : str
        Sample sheet text.
    *args : str
        Extra Snakemake arguments.
    config : str, optional
        YAML merged over the fixture config.

    Returns
    -------
    subprocess.CompletedProcess
        The finished Snakemake process.
    """
    (tmp_path / "samples.tsv").write_text(sheet)
    (tmp_path / "extra.yaml").write_text(f"samples: samples.tsv\n{config}")
    (tmp_path / "resources").symlink_to(ROOT / "resources")
    (tmp_path / ".tests").symlink_to(ROOT / ".tests")
    return snakemake(
        "--snakefile",
        ROOT / "workflow/Snakefile",
        "--directory",
        tmp_path,
        "--configfile",
        CONFIG,
        "extra.yaml",
        "-n",
        *args,
        cwd=tmp_path,
    )


@pytest.mark.parametrize(
    ("sheet", "message"),
    [
        ("experiment\treplicate\tr1\n", r"missing columns \['type'\]"),
        (HEADER, "no rows"),
        (
            HEADER + CHIP.replace("tf", "chip") + INPUT,
            "type must be tf, histone or control",
        ),
        (HEADER + CHIP.replace("Input", ""), "tf experiments need a control"),
        (
            HEADER + CHIP + INPUT.replace("control", "histone", 1),
            "must be an experiment of type control",
        ),
        (
            HEADER + CHIP + INPUT.replace("\t\t\t", "\tChip\t\t"),
            "cannot have a control",
        ),
        (HEADER + CHIP.replace("\t1\t", "\t2\t") + INPUT, "numbered 1 to N"),
        (HEADER + CHIP.replace("\t1\t", "\tone\t") + INPUT, "positive integer"),
        (HEADER + CHIP.replace("Chip", "Chip.1") + INPUT, "must match"),
        (HEADER + CHIP.replace(".gz", "") + INPUT, "must be gzipped"),
        (
            HEADER + CHIP + CHIP.replace("\t\n", "\tchip_R2.fastq.gz\n") + INPUT,
            "mix single and paired end",
        ),
        (
            HEADER + CHIP + CHIP.replace("tf\tInput", "histone\tInput") + INPUT,
            "rows disagree",
        ),
        (
            HEADER
            + CHIP
            + INPUT
            + INPUT.replace("\t1\t", "\t2\t").replace("\t\n", "\tinput_R2.fastq.gz\n"),
            "all replicates single or all paired end",
        ),
    ],
)
def test_sheet_errors(tmp_path, sheet, message):
    """Check that a bad sheet stops the run with a message naming the problem."""
    done = dry_run(tmp_path, sheet)
    assert done.returncode != 0
    assert re.search(message, done.stderr + done.stdout), done.stderr[-2000:]


def test_odd_paired_end_subsample_is_rejected(tmp_path):
    """Check that a paired-end library cannot be subsampled to an odd read count."""
    sheet = HEADER + CHIP.replace("\t\n", "\tchip_R2.fastq.gz\n") + INPUT
    done = dry_run(tmp_path, sheet, config="params:\n  subsample_reads: 101\n")
    assert done.returncode != 0
    assert "even number of reads" in done.stderr + done.stdout


def test_valid_sheet_plans_a_run(tmp_path):
    """Check that a minimal tf experiment with a control resolves to a job graph."""
    for fastq in ("chip_R1.fastq.gz", "input_R1.fastq.gz"):
        (tmp_path / fastq).touch()
    done = dry_run(tmp_path, HEADER + CHIP + INPUT)
    assert done.returncode == 0, done.stderr[-3000:]
    assert "choose_ctl" in done.stdout


FIXTURE_SHEET = ROOT / ".tests/integration/samples.tsv"
PREFIX = "results-fixture/fixture"


def test_tagalign_job_fails_on_a_corrupt_bam(tmp_path):
    """Check that a corrupt BAM fails the tagAlign job instead of leaving empty files."""
    done = run_rule(
        tmp_path,
        "tagalign",
        f"{PREFIX}/YCR/tagalign/YCR.rep1.tagAlign.gz",
        {f"{PREFIX}/YCR/filter/YCR.rep1.bam": b"not a BAM"},
    )
    assert done.returncode != 0


def test_frip_job_fails_on_unreadable_peaks(tmp_path, run):
    """Check that a bedtools error fails the FRiP job instead of reporting 0."""
    files = {
        f"{PREFIX}/BGT1/{p}": run / f"BGT1/{p}"
        for p in (
            "tagalign/BGT1.rep1.tagAlign.gz",
            "tagalign/BGT1.rep1.tagAlign.counts",
            "xcor/BGT1.rep1.fraglen",
            "overlap/BGT1.rep1-pr1_vs_rep1-pr2.bfilt.narrowPeak.gz",
        )
    }
    files["results-fixture/refs/reference.chromsizes"] = CHRSZ
    files[f"{PREFIX}/BGT1/macs3/BGT1.rep1.bfilt.narrowPeak.gz"] = gzip.compress(
        b"chr19\tnot\ta peak\n"
    )
    done = run_rule(tmp_path, "frip", f"{PREFIX}/BGT1/frip/BGT1.rep1.tsv", files)
    assert done.returncode != 0


def test_idr_with_too_few_peaks_writes_empty_sets(tmp_path, run):
    """Check that IDR on fewer than 20 merged peaks gives empty results, not a failure."""
    files = {"results-fixture/refs/reference.chromsizes": CHRSZ}
    for sample in ("rep1", "rep2", "pooled"):
        peaks = lines(run / f"YCR/spp/YCR.{sample}.regionPeak.gz")[:5]
        files[f"{PREFIX}/YCR/spp/YCR.{sample}.regionPeak.gz"] = gzip.compress(
            "".join(f"{x}\n" for x in peaks).encode()
        )
    target = f"{PREFIX}/YCR/idr/YCR.rep1_vs_rep2.regionPeak.gz"
    done = run_rule(tmp_path, "idr", target, files)
    assert done.returncode == 0, done.stderr[-3000:]
    assert lines(tmp_path / target) == []


def test_shared_control_is_processed_once(run):
    """Check that two experiments sharing a control call against the same files."""
    tables = [
        pl.read_csv(run / f"{exp}/choose_ctl/{exp}.tsv", separator="\t")
        for exp in ("YCR", "YCR_h")
    ]
    assert tables[0].equals(tables[1])
    subsamples = sorted((run / "YCR_ctl/tagalign").glob("*.sub*.tagAlign.gz"))
    assert len(subsamples) == tables[0]["subsample"].n_unique()


def test_xcor_trims_only_when_every_read_is_long_enough(run):
    """Check ENCODE's all-or-nothing trim: 100 bp YCR reads trimmed, 36 bp BGT reads not."""
    note = "no trimming"
    assert note not in (ROOT / "logs/fixture/xcor_tagalign/YCR.rep1.txt").read_text()
    assert note in (ROOT / "logs/fixture/xcor_tagalign/BGT.rep1.txt").read_text()


def test_true_rep_only_drops_pseudoreplicates(tmp_path):
    """Check that true_rep_only plans no pseudoreplicate or reproducibility jobs."""
    done = dry_run(
        tmp_path,
        FIXTURE_SHEET.read_text(),
        "--forceall",
        config="params:\n  true_rep_only: true\n",
    )
    assert done.returncode == 0, done.stderr[-3000:]
    jobs = done.stdout
    for rule in ("pseudoreps", "macs3_pseudorep", "reproducibility"):
        assert f"\n{rule} " not in jobs
    assert "\noverlap " in jobs and "\nidr " in jobs
