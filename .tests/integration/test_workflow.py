"""Check sheet validation, failure behavior and the shape of the job graph."""

import re

import pytest
from conftest import CONFIG, ROOT, snakemake

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
