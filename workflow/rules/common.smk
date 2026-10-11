# Config, experiment sheet, and shared helpers.

import re
from itertools import combinations
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
TMP_DIR = config["tmp_dir"]

REFERENCE_FASTA = config["references"]["fasta"]
BLACKLIST = config["references"]["blacklist"]
CHROMSIZES = f"{REFDIR}/reference.chromsizes"
INDEX = {
    "bowtie2": f"{REFDIR}/bowtie2_build/reference",
    "bwa": f"{REFDIR}/bwa_index/reference",
}
INDEX_FILES = {
    "bowtie2": multiext(
        INDEX["bowtie2"], ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"
    ),
    "bwa": multiext(INDEX["bwa"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
}
ALIGNER = config["params"]["aligner"]
if ALIGNER not in INDEX:
    raise WorkflowError(f"params.aligner must be one of {sorted(INDEX)}, not {ALIGNER!r}.")

# The workspace path is relative to the rule file. Every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)
IDR_ENV = pixi(workspace="../envs", env="idr", locked=True)
SPP_ENV = pixi(workspace="../envs", env="spp", locked=True)

# samtools sort -T and the per-job scratch directories live here.
Path(TMP_DIR).mkdir(parents=True, exist_ok=True)

# Dots separate name parts in file names, so names never contain one.
NAME = r"[A-Za-z0-9][A-Za-z0-9_-]*"
CALLER = {"tf": "spp", "histone": "macs3"}
PEAK_TYPE = {"spp": "regionPeak", "macs3": "narrowPeak"}


_SHEET = Path(config["samples"])
if not _SHEET.exists():
    raise WorkflowError(
        f"Sample sheet not found: {_SHEET}\n"
        "Columns: experiment, replicate, type, r1, [r2], [control], [crop_length]; "
        "see config/samples.example.tsv."
    )


def _fail(message):
    """Raise a sample sheet error.

    Parameters
    ----------
    message : str
        What is wrong.

    Raises
    ------
    WorkflowError
        Always, prefixed with the sheet path.
    """
    raise WorkflowError(f"{_SHEET}: {message}")


# All-string read; optional columns default to "".
_sheet = pl.read_csv(_SHEET, separator="\t", infer_schema_length=0, comment_prefix="#")
_missing = {"experiment", "replicate", "type", "r1"} - set(_sheet.columns)
if _missing:
    _fail(f"missing columns {sorted(_missing)}")
_sheet = _sheet.with_columns(
    pl.lit("").alias(c) for c in ("r2", "control", "crop_length") if c not in _sheet.columns
).with_columns(pl.all().fill_null("").str.strip_chars())
if _sheet.is_empty():
    _fail("no rows")

EXPERIMENTS = {}  # exp -> {"type", "control", "crop", "libs"}
LIBRARIES = {}  # (exp, lib) -> {"r1", "r2", "layout"}
for row in _sheet.iter_rows(named=True):
    exp = row["experiment"]
    where = f"{exp} replicate {row['replicate']}"
    if not re.fullmatch(NAME, exp):
        _fail(f"experiment {exp!r} must match {NAME}")
    if not re.fullmatch(r"[1-9][0-9]*", row["replicate"]):
        _fail(f"{where}: replicate must be a positive integer")
    if row["type"] not in ("tf", "histone", "control"):
        _fail(f"{where}: type must be tf, histone or control, not {row['type']!r}")
    if not re.fullmatch(r"[0-9]*", row["crop_length"]):
        _fail(f"{where}: crop_length must be a non-negative integer")
    if not row["r1"]:
        _fail(f"{where}: empty r1")
    for fastq in (row["r1"], row["r2"]):
        if fastq and not fastq.endswith(".gz"):
            _fail(f"{where}: {fastq} must be gzipped (.gz)")
    attributes = {
        "type": row["type"],
        "control": row["control"] or None,
        "crop": int(row["crop_length"] or 0),
    }
    e = EXPERIMENTS.setdefault(exp, {**attributes, "libs": []})
    if {k: e[k] for k in attributes} != attributes:
        _fail(f"{exp}: rows disagree on type, control or crop_length")
    lib = f"{'ctl' if row['type'] == 'control' else 'rep'}{int(row['replicate'])}"
    layout = "pe" if row["r2"] else "se"
    entry = LIBRARIES.get((exp, lib))
    if entry is None:
        e["libs"].append(lib)
        entry = LIBRARIES[exp, lib] = {"r1": [], "r2": [], "layout": layout}
    if entry["layout"] != layout:
        _fail(f"{where}: lanes mix single and paired end")
    entry["r1"].append(row["r1"])
    if row["r2"]:
        entry["r2"].append(row["r2"])

# Control subsampling treats a pooled control as having the layout of its first replicate.
_CTL_SUBSAMPLING = (
    config["params"]["ctl_depth_limit"] > 0
    or config["params"]["exp_ctl_depth_ratio_limit"] > 0
    or config["params"]["ctl_subsample_reads"] > 0
)
for exp, e in EXPERIMENTS.items():
    e["libs"].sort(key=lambda lib: int(lib[3:]))
    if [int(lib[3:]) for lib in e["libs"]] != list(range(1, len(e["libs"]) + 1)):
        _fail(f"{exp}: replicates must be numbered 1 to N")
    layouts = {LIBRARIES[exp, lib]["layout"] for lib in e["libs"]}
    if e["type"] == "control":
        if e["control"]:
            _fail(f"{exp}: a control experiment cannot have a control")
        if _CTL_SUBSAMPLING and len(layouts) > 1:
            _fail(f"{exp}: control subsampling needs all replicates single or all paired end")
        n = config["params"]["ctl_subsample_reads"]
    else:
        if e["type"] == "tf" and not e["control"]:
            _fail(f"{exp}: tf experiments need a control, which SPP calls against")
        ctl = EXPERIMENTS.get(e["control"]) if e["control"] else None
        if e["control"] and (ctl is None or ctl["type"] != "control"):
            _fail(f"{exp}: control {e['control']!r} must be an experiment of type control")
        n = config["params"]["subsample_reads"]
    if n % 2 and "pe" in layouts:
        _fail(f"{exp}: paired-end subsampling needs an even number of reads, not {n}")

CHIPS = [exp for exp, e in EXPERIMENTS.items() if e["type"] != "control"]


wildcard_constraints:
    exp=NAME,
    lib=r"(rep|ctl)[1-9][0-9]*",
    rep=r"rep[1-9][0-9]*",
    base=r"rep[1-9][0-9]*|pooled",
    pr=r"pr[12]",
    sample=r"(rep[1-9][0-9]*|pooled)(-pr[12])?",
    pool=r"pooled(-pr[12])?",
    name=r"((rep|ctl)[1-9][0-9]*|pooled)(-pr[12])?(\.sub[1-9][0-9]*)*",
    n=r"[1-9][0-9]*",
    lane=r"[1-9][0-9]*",
    pair=r"[A-Za-z0-9-]+_vs_[A-Za-z0-9-]+",
    method=r"idr|overlap",
    caller=r"spp|macs3",
    peakdir=r"spp|macs3|idr|overlap",
    unit=r"[A-Za-z0-9_-]+",
    ptype=r"narrowPeak|regionPeak",
    track=r"fc|pval",
    strand=r"positive|negative",


# Shell functions shared by the tagAlign rules. count passes lines through and writes
# their number and byte count to $1, so no rule reads a tagAlign twice to size it.
# keystream is ENCODE's seeded shuffle source; the ENCODE image's OpenSSL 1.0.2
# derives the key with MD5, which -md md5 keeps on newer OpenSSL.
TA_FUNCTIONS = r"""
        count() {{ awk -v out="$1" '{{ print; n++; b += length($0) + 1 }} END {{ printf "%d\t%d\n", n, b > out }}'; }}
        keystream() {{ openssl enc -aes-256-ctr -md md5 -pass pass:"$1" -nosalt < /dev/zero 2> /dev/null; }}
        unpair() {{ awk 'BEGIN {{ OFS = "\t" }} {{ print $1, $2, $3, $4, $5, $6; print $7, $8, $9, $10, $11, $12 }}'; }}
"""


def lanes(exp, lib, mate):
    """Return a library's FASTQs of one mate, in sheet order.

    Parameters
    ----------
    exp : str
        Experiment.
    lib : str
        Library, `repN` or `ctlN`.
    mate : {"r1", "r2"}
        Mate.

    Returns
    -------
    list of str
        FASTQ paths. Empty for `r2` of a single-end library.
    """
    return LIBRARIES[exp, lib][mate]


def layout(exp, name):
    """Return the layout of a library or of the tagAlign built from it.

    Parameters
    ----------
    exp : str
        Experiment.
    name : str
        tagAlign name. A pool takes the layout of the experiment's first library.

    Returns
    -------
    {"se", "pe"}
        Single or paired end.
    """
    lib = re.match(r"(rep|ctl)[0-9]+", name)
    return LIBRARIES[exp, lib[0] if lib else EXPERIMENTS[exp]["libs"][0]]["layout"]


def ta(exp, name):
    """Return the path of a tagAlign.

    Parameters
    ----------
    exp : str
        Experiment.
    name : str
        tagAlign name, such as `rep1`, `rep1.sub20000000`, `pooled` or `rep1-pr2`.

    Returns
    -------
    str
        Path `<outdir>/<run_id>/<exp>/tagalign/<exp>.<name>.tagAlign.gz`.
    """
    return f"{OUTDIR}/{exp}/tagalign/{exp}.{name}.tagAlign.gz"


def counts(exp, name):
    """Return the path of a tagAlign's line and byte counts.

    Parameters
    ----------
    exp : str
        Experiment.
    name : str
        tagAlign name.

    Returns
    -------
    str
        Path of the one-line `lines<TAB>bytes` file beside the tagAlign.
    """
    return f"{OUTDIR}/{exp}/tagalign/{exp}.{name}.tagAlign.counts"


def effective(exp, lib):
    """Return the tagAlign name a library contributes downstream.

    Parameters
    ----------
    exp : str
        Experiment.
    lib : str
        Library, `repN` or `ctlN`.

    Returns
    -------
    str
        `lib`, or `lib.sub<N>` when `subsample_reads` (ChIP) or
        `ctl_subsample_reads` (control) is set.
    """
    key = "ctl_subsample_reads" if lib.startswith("ctl") else "subsample_reads"
    n = config["params"][key]
    return f"{lib}.sub{n}" if n else lib


def pool_members(exp, pool):
    """Return the tagAlign names a pool concatenates.

    Parameters
    ----------
    exp : str
        Experiment.
    pool : str
        `pooled`, `pooled-pr1` or `pooled-pr2`.

    Returns
    -------
    list of str
        One tagAlign name per library.
    """
    suffix = pool.removeprefix("pooled")
    return [lib + suffix if suffix else effective(exp, lib) for lib in EXPERIMENTS[exp]["libs"]]


def bases(exp):
    """Return the replicates of an experiment, plus `pooled` when there are two or more.

    Parameters
    ----------
    exp : str
        ChIP experiment.

    Returns
    -------
    list of str
        `rep1` to `repN`, then `pooled` if N >= 2.
    """
    reps = EXPERIMENTS[exp]["libs"]
    return reps + ["pooled"] * (len(reps) > 1)


def samples(exp):
    """Return every tagAlign an experiment calls peaks on.

    Parameters
    ----------
    exp : str
        ChIP experiment.

    Returns
    -------
    list of str
        `bases(exp)`, then each one's pseudoreplicates unless `true_rep_only`.
    """
    if config["params"]["true_rep_only"]:
        return bases(exp)
    return bases(exp) + [f"{b}-pr{k}" for b in bases(exp) for k in (1, 2)]


def sample_ta(exp, sample):
    """Return the tagAlign name a peak-calling sample reads.

    Parameters
    ----------
    exp : str
        ChIP experiment.
    sample : str
        `repN`, `pooled`, or a pseudoreplicate of either.

    Returns
    -------
    str
        `effective(exp, sample)` for a replicate, else `sample`.
    """
    return effective(exp, sample) if sample.startswith("rep") and "-" not in sample else sample


def pairs(exp):
    """Return the peak set pairs compared by overlap and IDR, in ENCODE's order.

    Parameters
    ----------
    exp : str
        ChIP experiment.

    Returns
    -------
    list of str
        True replicate pairs `repI_vs_repJ`, then unless `true_rep_only` each
        replicate's `repN-pr1_vs_repN-pr2`, then `pooled-pr1_vs_pooled-pr2`
        when there are two or more replicates.
    """
    reps = EXPERIMENTS[exp]["libs"]
    true = [f"{a}_vs_{b}" for a, b in combinations(reps, 2)]
    if config["params"]["true_rep_only"]:
        return true
    return true + [f"{b}-pr1_vs_{b}-pr2" for b in bases(exp)]


def oracle(pair):
    """Return the sample whose peaks anchor a pair, as ENCODE's chip.wdl wires it.

    Its tagAlign and fragment length also score the pair's FRiP.

    Parameters
    ----------
    pair : str
        Pair from `pairs()`.

    Returns
    -------
    str
        `repN` for `repN-pr1_vs_repN-pr2`, else `pooled`.
    """
    first = pair.split("_vs_")[0]
    return first.split("-")[0] if "-" in first else "pooled"


def methods(exp):
    """Return the reproducibility methods of an experiment.

    Parameters
    ----------
    exp : str
        ChIP experiment.

    Returns
    -------
    list of str
        `overlap`, plus `idr` for tf experiments.
    """
    return ["idr", "overlap"] if EXPERIMENTS[exp]["type"] == "tf" else ["overlap"]


def peak_type(exp):
    """Return the peak file type of an experiment.

    Parameters
    ----------
    exp : str
        ChIP experiment.

    Returns
    -------
    str
        `regionPeak` for SPP (tf), `narrowPeak` for MACS3 (histone).
    """
    return PEAK_TYPE[CALLER[EXPERIMENTS[exp]["type"]]]


def peaks(exp, sample, bfilt=False):
    """Return the path of a sample's peaks from the experiment's peak caller.

    Parameters
    ----------
    exp : str
        ChIP experiment.
    sample : str
        Peak-calling sample.
    bfilt : bool, default False
        Return the exclusion-list-filtered peaks.

    Returns
    -------
    str
        Path under `spp/` or `macs3/`.
    """
    caller = CALLER[EXPERIMENTS[exp]["type"]]
    return f"{OUTDIR}/{exp}/{caller}/{exp}.{sample}.{'bfilt.' * bfilt}{peak_type(exp)}.gz"


def pair_peaks(exp, method, pair, bfilt=True):
    """Return the path of a pair's overlap or IDR peaks.

    Parameters
    ----------
    exp : str
        ChIP experiment.
    method : {"idr", "overlap"}
        Reproducibility method.
    pair : str
        Pair from `pairs()`.
    bfilt : bool, default True
        Return the exclusion-list-filtered peaks.

    Returns
    -------
    str
        Path under `idr/` or `overlap/`.
    """
    return f"{OUTDIR}/{exp}/{method}/{exp}.{pair}.{'bfilt.' * bfilt}{peak_type(exp)}.gz"


def fraglen(exp, sample):
    """Return the fragment length file a sample's peak call and FRiP use.

    Parameters
    ----------
    exp : str
        ChIP experiment.
    sample : str
        Peak-calling sample. Pseudoreplicates use their replicate's length,
        pools the rounded mean.

    Returns
    -------
    str
        Path under `xcor/`.
    """
    return f"{OUTDIR}/{exp}/xcor/{exp}.{sample.split('-')[0]}.fraglen"


def control_ta(exp, sample):
    """Return the control tagAlign a sample calls peaks against.

    Reads the experiment's `choose_ctl` checkpoint table.

    Parameters
    ----------
    exp : str
        ChIP experiment.
    sample : str
        Peak-calling sample. Pseudoreplicates use their replicate's control.

    Returns
    -------
    list of str
        One tagAlign path, or none for an experiment without a control.
    """
    ctl_exp = EXPERIMENTS[exp]["control"]
    if ctl_exp is None:
        return []
    table = pl.read_csv(checkpoints.choose_ctl.get(exp=exp).output[0], separator="\t")
    row = table.row(by_predicate=pl.col("sample") == sample.split("-")[0], named=True)
    name = row["control"] if row["control"] == "pooled" else effective(ctl_exp, row["control"])
    return [ta(ctl_exp, f"{name}.sub{row['subsample']}" if row["subsample"] else name)]


def frip_sets(exp, sample):
    """Return the filtered peak sets whose FRiP a sample's tagAlign scores.

    Parameters
    ----------
    exp : str
        ChIP experiment.
    sample : str
        Peak-calling sample.

    Returns
    -------
    dict of str to str
        Path by label: the caller's name for the sample's own peaks, and
        `<method>.<pair>` for each pair the sample anchors.
    """
    sets = {CALLER[EXPERIMENTS[exp]["type"]]: peaks(exp, sample, bfilt=True)}
    for pair in pairs(exp):
        if oracle(pair) == sample:
            for method in methods(exp):
                sets[f"{method}.{pair}"] = pair_peaks(exp, method, pair)
    return sets
