# Config, sample sheet, and shared helpers.

from pathlib import Path

import polars as pl
from snakemake.exceptions import WorkflowError


configfile: "config/config.yaml"


RESULTS = config["outdir"]
RUN_ID = config.get("run_id", "default")
OUTDIR = f"{RESULTS}/{RUN_ID}"          # per-run
REFDIR = f"{RESULTS}/refs"              # run-independent
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"
REFERENCE = config["references"]["reference"]

# The workspace path is relative to the rule file. Every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)


_SHEET = Path(config["samples"])
if not _SHEET.exists():
    raise WorkflowError(
        f"Sample sheet not found: {_SHEET}\n"
        "Columns: sample_id, input, [label]; see config/samples.example.tsv."
    )

# All-string read so empty cells are "" (absent), not a type error.
_manifest = pl.read_csv(_SHEET, separator="\t", infer_schema_length=0, comment_prefix="#")


def _column(name):
    """Map each sample to its value in sheet column `name`.

    Parameters
    ----------
    name : str
        Sample sheet column.

    Returns
    -------
    dict of str to str or None
        Value per `sample_id`, with empty cells as None. Empty when the sheet
        has no such column.
    """
    if name not in _manifest.columns:
        return {}
    return {
        row["sample_id"]: (row[name] or None)
        for row in _manifest.iter_rows(named=True)
    }


SAMPLES = _manifest.get_column("sample_id").to_list()
INPUT_OF = _column("input")
LABEL_OF = _column("label")

if not SAMPLES:
    raise WorkflowError(f"Sample sheet {_SHEET} has no rows.")


wildcard_constraints:
    sample=r"[A-Za-z0-9][A-Za-z0-9_.-]*",


def prefix(sample):
    """Return the output prefix of a sample.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.

    Returns
    -------
    str
        Path `<outdir>/<run_id>/<sample>/<sample>`.
    """
    return f"{OUTDIR}/{sample}/{sample}"


def analyze_flags():
    """Return the analyze.py flag tokens from `config["analyze"]`.

    Returns
    -------
    list of str
        Flag tokens. Interpolate them with `:q`.
    """
    a = config["analyze"]
    flags = [
        "--window", str(a["window"]),
        "--threshold", str(a["threshold"]),
        "--metric", a["metric"],
    ]
    if a["normalize"]:
        flags.append("--normalize")
    if a["random_state"] is not None:
        flags += ["--random_state", str(a["random_state"])]
    return flags
