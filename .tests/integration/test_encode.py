"""Compare each workflow step with ENCODE's own task script on the same inputs."""

import subprocess

import polars as pl
import pytest
from conftest import (
    BLACKLIST,
    CHRSZ,
    FIXTURE,
    GENOME,
    REFS,
    ROOT,
    bigwig_stats,
    encode,
    first_difference,
    lines,
)

# Library, layout and ENCODE's paired-end flag of the two fixture experiments.
LIBRARIES = [("YCR", "rep1", "--paired-end"), ("BGT", "rep1", "")]
REGEX = "--regex-bfilt-peak-chr-name 'chr[\\dXY]+'"


def alignments(bam):
    """List a BAM's records by name, flag, position, MAPQ and CIGAR.

    Parameters
    ----------
    bam : pathlib.Path
        BAM file.

    Returns
    -------
    list of str
        Sorted first six SAM fields of every record.
    """
    sam = subprocess.run(
        ["samtools", "view", str(bam)], capture_output=True, text=True, check=True
    ).stdout
    return sorted("\t".join(r.split("\t")[:6]) for r in sam.splitlines())


def jaccard(a, b):
    """Return the Jaccard index of two collections as sets.

    Parameters
    ----------
    a, b : iterable
        Collections to compare.

    Returns
    -------
    float
        Shared over combined distinct items, 1 when both are empty.
    """
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else 1.0


def regions(path):
    """Return the regions (contig, start, end) of a peak file.

    Parameters
    ----------
    path : pathlib.Path
        Peak file, gzipped or not.

    Returns
    -------
    set of tuple of str
        One tuple per line.
    """
    return {tuple(line.split("\t")[:3]) for line in lines(path)}


def frip(run, exp, sample, peak_set):
    """Return one FRiP row of a sample's table as reads in peaks over reads.

    Parameters
    ----------
    run : pathlib.Path
        Workflow run directory.
    exp, sample : str
        Experiment and peak-calling sample.
    peak_set : str
        Row label: the caller, or ``<method>.<pair>``.

    Returns
    -------
    float
        Fraction of tagAlign lines in the peak set.
    """
    table = pl.read_csv(run / f"{exp}/frip/{exp}.{sample}.tsv", separator="\t")
    row = table.row(by_predicate=pl.col("peak_set") == peak_set, named=True)
    return row["reads_in_peaks"] / row["reads"]


def control_of(run, exp, sample):
    """Return the control tagAlign the workflow called a sample against.

    Parameters
    ----------
    run : pathlib.Path
        Workflow run directory.
    exp : str
        ChIP experiment.
    sample : str
        Replicate or ``pooled``.

    Returns
    -------
    pathlib.Path
        Control tagAlign, from the experiment's control choice table.
    """
    table = pl.read_csv(run / f"{exp}/choose_ctl/{exp}.tsv", separator="\t")
    row = table.row(by_predicate=pl.col("sample") == sample, named=True)
    ctl_exp = f"{exp.split('_')[0]}_ctl"
    name = row["control"] + (f".sub{row['subsample']}" if row["subsample"] else "")
    return run / f"{ctl_exp}/tagalign/{ctl_exp}.{name}.tagAlign.gz"


def fraglen(run, exp, sample):
    """Return the fragment length the workflow used for a sample.

    Parameters
    ----------
    run : pathlib.Path
        Workflow run directory.
    exp, sample : str
        Experiment and replicate or ``pooled``.

    Returns
    -------
    int
        Fragment length.
    """
    return int(lines(run / f"{exp}/xcor/{exp}.{sample}.fraglen")[0])


def encode_filter(run, exp, lib, paired):
    """Filter a raw BAM with ENCODE's unique-read, sambamba path.

    Parameters
    ----------
    run : pathlib.Path
        Workflow run directory.
    exp, lib : str
        Experiment and library.
    paired : str
        ``--paired-end`` or empty.

    Returns
    -------
    pathlib.Path
        ENCODE output directory.
    """
    bam = run / f"{exp}/align/{exp}.{lib}.bam"
    # Its last step copies the input BAM into --out-dir, which must not hold the link.
    out = encode(
        f"filter/{exp}.{lib}",
        [bam, CHRSZ],
        f"python3 $(which encode_task_filter.py) {bam.name} {paired} "
        "--multimapping 0 --dup-marker sambamba --mapq-thresh 30 --filter-chrs "
        f"--chrsz {CHRSZ.name} --mito-chr-name chrM --nth 8 --mem-gb 8 --out-dir out",
    )
    return out / "out"


@pytest.mark.parametrize(("exp", "lib", "paired"), LIBRARIES)
def test_raw_alignment_counts_match_encode(run, exp, lib, paired):
    """Check total and mapped reads of the raw BAM against ENCODE's SAMstats."""
    bam = run / f"{exp}/align/{exp}.{lib}.bam"
    fastq = FIXTURE / f"{exp}_{lib}_R1.fastq.gz"
    out = encode(
        f"post_align/{exp}.{lib}",
        [bam, fastq],
        f"python3 $(which encode_task_post_align.py) {fastq.name} {bam.name} "
        "--mito-chr-name chrM --nth 4 --mem-gb 4",
    )
    samstats = lines(next(out.glob("*.samstats.qc")))
    ours = pl.read_json(run / f"{exp}/align/{exp}.{lib}.flagstat.json")
    passed = ours["QC-passed reads"][0]
    assert int(samstats[0].split()[0]) == passed["total"]
    mapped = next(x for x in samstats if " mapped (" in x)
    assert int(mapped.split()[0]) == passed["mapped"]


@pytest.mark.parametrize(("exp", "lib", "paired"), LIBRARIES)
def test_filtered_bam_matches_encode(run, exp, lib, paired):
    """Check that filtering and duplicate removal keep ENCODE's records."""
    ours = alignments(run / f"{exp}/filter/{exp}.{lib}.bam")
    theirs = alignments(next(encode_filter(run, exp, lib, paired).glob("*.nodup.bam")))
    assert first_difference(ours, theirs) is None


@pytest.mark.parametrize(("exp", "lib", "paired"), LIBRARIES)
def test_duplicates_match_encode(run, exp, lib, paired):
    """Check the duplicate count against sambamba's log in ENCODE's run."""
    log = "\n".join(lines(encode_filter(run, exp, lib, paired).parent / "{}.dup.qc"))
    found = int(log.split("found ")[1].split()[0])
    ours = pl.read_json(run / f"{exp}/filter/{exp}.{lib}.dup.flagstat.json")
    assert ours["QC-passed reads"][0]["duplicates"] == found


@pytest.mark.parametrize(("exp", "lib", "paired"), LIBRARIES)
def test_library_complexity_matches_encode(run, exp, lib, paired):
    """Check NRF, PBC1 and PBC2 and their counts against ENCODE's."""
    qc = next(encode_filter(run, exp, lib, paired).glob("*.lib_complexity.qc"))
    mt, m0, m1, m2, nrf, pbc1, pbc2 = lines(qc)[0].split("\t")
    ours = pl.read_json(run / f"{exp}/pbc_qc/{exp}.{lib}.json").row(0, named=True)
    assert ours["total_fragments"] == int(mt)
    assert ours["distinct_fragments"] == int(m0)
    assert ours["positions_with_one_read"] == int(m1)
    assert ours["positions_with_two_reads"] == int(m2)
    for key, value in (("NRF", nrf), ("PBC1", pbc1), ("PBC2", pbc2)):
        assert ours[key] == pytest.approx(float(value), abs=1e-6)


def encode_tagalign(run, exp, lib, paired):
    """Convert a filtered BAM to a tagAlign with ENCODE's bam2ta.

    Parameters
    ----------
    run : pathlib.Path
        Workflow run directory.
    exp, lib : str
        Experiment and library.
    paired : str
        ``--paired-end`` or empty.

    Returns
    -------
    pathlib.Path
        ENCODE's tagAlign.
    """
    bam = run / f"{exp}/filter/{exp}.{lib}.bam"
    out = encode(
        f"bam2ta/{exp}.{lib}",
        [bam],
        f"python3 $(which encode_task_bam2ta.py) {bam.name} {paired} "
        "--disable-tn5-shift --mito-chr-name chrM --subsample 0 --nth 4",
    )
    return out / f"{exp}.{lib}.tagAlign.gz"


@pytest.mark.parametrize(("exp", "lib", "paired"), LIBRARIES)
def test_tagalign_matches_encode(run, exp, lib, paired):
    """Check that the tagAlign has ENCODE's lines in ENCODE's order, and its counts."""
    ours = lines(run / f"{exp}/tagalign/{exp}.{lib}.tagAlign.gz")
    assert first_difference(ours, lines(encode_tagalign(run, exp, lib, paired))) is None
    count, size = map(
        int, lines(run / f"{exp}/tagalign/{exp}.{lib}.tagAlign.counts")[0].split()
    )
    assert (count, size) == (len(ours), sum(len(x) + 1 for x in ours))


@pytest.mark.parametrize(("exp", "lib", "paired"), LIBRARIES)
def test_pseudoreplicates_match_encode(run, exp, lib, paired):
    """Check that pseudoreplicates split the reads (whole pairs) as ENCODE's spr."""
    ta = run / f"{exp}/tagalign/{exp}.{lib}.tagAlign.gz"
    out = encode(
        f"spr/{exp}.{lib}",
        [ta],
        f"python3 $(which encode_task_spr.py) {ta.name} {paired} "
        "--pseudoreplication-random-seed 0",
    )
    for pr in (1, 2):
        ours = lines(run / f"{exp}/tagalign/{exp}.{lib}-pr{pr}.tagAlign.gz")
        theirs = lines(out / f"{exp}.{lib}.pr{pr}.tagAlign.gz")
        assert first_difference(ours, theirs) is None


@pytest.mark.parametrize("exp", ["YCR", "BGT_ctl"])
def test_pool_matches_encode(run, exp):
    """Check a pooled tagAlign against ENCODE's pool_ta."""
    tas = sorted((run / f"{exp}/tagalign").glob(f"{exp}.[rc][et][pl][0-9].tagAlign.gz"))
    out = encode(
        f"pool_ta/{exp}",
        tas,
        f"python3 $(which encode_task_pool_ta.py) {' '.join(t.name for t in tas)} --prefix x",
    )
    ours = lines(run / f"{exp}/tagalign/{exp}.pooled.tagAlign.gz")
    assert first_difference(ours, lines(out / "x.pooled.tagAlign.gz")) is None


@pytest.mark.parametrize(
    ("exp", "ctl_exp", "paired"),
    [("YCR", "YCR_ctl", "--paired-end"), ("BGT", "BGT_ctl", "")],
)
def test_control_subsample_matches_encode(run, exp, ctl_exp, paired):
    """Check the subsampled pooled control against ENCODE's subsample_ctl."""
    control = control_of(run, exp, "pooled")
    assert ".sub" in control.name, "the fixture config should force a subsample"
    n = control.name.split(".sub")[1].split(".")[0]
    ta = run / f"{ctl_exp}/tagalign/{ctl_exp}.pooled.tagAlign.gz"
    out = encode(
        f"subsample_ctl/{ctl_exp}.{n}",
        [ta],
        f"python3 $(which encode_task_subsample_ctl.py) {ta.name} {paired} --subsample {n}",
    )
    theirs = lines(next(out.glob("*.pooled.*.tagAlign.gz")))
    assert first_difference(lines(control), theirs) is None


@pytest.mark.parametrize(
    "flags",
    [
        "--always-use-pooled-ctl --ctl-depth-ratio 1.2 --ctl-depth-limit 200000000 --exp-ctl-depth-ratio-limit 5.0",
        "--always-use-pooled-ctl --ctl-depth-ratio 1.2 --ctl-depth-limit 46000 --exp-ctl-depth-ratio-limit 0.9",
        "--ctl-depth-ratio 1.2 --ctl-depth-limit 46000 --exp-ctl-depth-ratio-limit 0.9",
        "--ctl-depth-ratio 100 --ctl-depth-limit 0 --exp-ctl-depth-ratio-limit 0",
    ],
)
@pytest.mark.parametrize(
    ("exp", "ctl_exp", "pe"), [("YCR", "YCR_ctl", True), ("BGT", "BGT_ctl", False)]
)
def test_control_choice_matches_encode(run, tmp_path, exp, ctl_exp, pe, flags):
    """Check control choice and subsample depths against ENCODE's choose_ctl."""
    reps = sorted((run / f"{exp}/tagalign").glob(f"{exp}.rep[0-9].tagAlign.gz"))
    ctls = sorted((run / f"{ctl_exp}/tagalign").glob(f"{ctl_exp}.ctl[0-9].tagAlign.gz"))
    out = encode(
        f"choose_ctl/{exp}.{abs(hash(flags))}",
        reps + ctls,
        f"python3 $(which encode_task_choose_ctl.py) --tas {' '.join(t.name for t in reps)} "
        f"--ctl-tas {' '.join(t.name for t in ctls)} {flags}",
    )
    table = tmp_path / "chosen.tsv"
    subprocess.run(
        [
            "python",
            ROOT / "workflow/scripts/choose_ctl.py",
            "--reps",
            *[str(t).replace(".gz", ".counts") for t in reps],
            "--controls",
            *[str(t).replace(".gz", ".counts") for t in ctls],
            "--output",
            table,
            *flags.split(),
            *(["--paired-end-control"] if pe else []),
        ],
        check=True,
    )
    ours = pl.read_csv(table, separator="\t")
    ids = [int(x) for x in lines(out / "chosen_ctl.tsv")]
    expected = ["pooled" if i < 0 else f"ctl{i + 1}" for i in ids]
    assert ours["control"].to_list()[:-1] == expected
    subsample = [int(x) for x in lines(out / "chosen_ctl_subsample.tsv")]
    subsample.append(int(lines(out / "chosen_ctl_subsample_pooled.txt")[0]))
    # ENCODE's subsample_ctl raises on an odd paired-end depth; ours keeps whole pairs.
    assert ours["subsample"].to_list() == [n - n % 2 if pe else n for n in subsample]


def encode_xcor_tagalign(run, exp, lib, paired):
    """Build xcor reads with ENCODE's chain: R1 trimmed, aligned, filtered, converted.

    Parameters
    ----------
    run : pathlib.Path
        Workflow run directory.
    exp, lib : str
        Experiment and library.
    paired : str
        Unused; xcor reads are single end.

    Returns
    -------
    pathlib.Path
        ENCODE's xcor tagAlign.
    """
    fastq = FIXTURE / f"{exp}_{lib}_R1.fastq.gz"
    index = sorted((REFS / "bowtie2_build").glob("reference.*.bt2l"))
    out = encode(
        f"xcor_reads/{exp}.{lib}",
        [fastq, *index],
        f"python3 $(which encode_task_trim_fastq.py) {fastq.name} --trim-bp 50 --out-dir trim && "
        "python3 $(which encode_task_bowtie2.py) reference trim/*.fastq.gz --nth 8 --mem-gb 8 --out-dir align && "
        "python3 $(which encode_task_filter.py) align/*.bam --no-dup-removal --multimapping 0 "
        "--dup-marker sambamba --mapq-thresh 30 --filter-chrs --mito-chr-name chrM --nth 8 --mem-gb 8 "
        "--out-dir filter && "
        "python3 $(which encode_task_bam2ta.py) filter/*.filt.bam --disable-tn5-shift --subsample 0 "
        "--mito-chr-name chrM --nth 4 --out-dir ta",
    )
    return next((out / "ta").glob("*.tagAlign.gz"))


@pytest.mark.parametrize(("exp", "lib", "paired"), LIBRARIES)
def test_xcor_reads_match_encode(run, exp, lib, paired):
    """Check the xcor tagAlign against ENCODE's trim, align and filter chain.

    bowtie2 2.5.5 and ENCODE's 2.3.4.3 may place a few reads differently.
    """
    ours = lines(run / f"{exp}/xcor/{exp}.{lib}.tagAlign.gz")
    theirs = lines(encode_xcor_tagalign(run, exp, lib, paired))
    assert jaccard(ours, theirs) >= 0.999


@pytest.mark.parametrize(
    ("exp", "kind"), [("YCR", "tf"), ("YCR_h", "histone"), ("BGT", "tf")]
)
def test_cross_correlation_matches_encode(run, exp, kind):
    """Check fragment length, NSC and RSC against ENCODE's xcor on the same reads."""
    ta = run / f"{exp}/xcor/{exp}.rep1.tagAlign.gz"
    out = encode(
        f"xcor/{exp}.rep1",
        [ta],
        f"python3 $(which encode_task_xcor.py) {ta.name} --mito-chr-name chrM "
        f"--subsample 15000000 --chip-seq-type {kind} --exclusion-range-min -500 --nth 2",
    )
    theirs = lines(next(out.glob("*.cc.qc")))[0].split("\t")
    ours = lines(run / f"{exp}/xcor/{exp}.rep1.cc.qc")[0].split("\t")
    assert ours[1:3] == theirs[1:3]
    assert float(ours[8]) == pytest.approx(float(theirs[8]), rel=1e-3)
    assert float(ours[9]) == pytest.approx(float(theirs[9]), rel=1e-3)


def test_pooled_fragment_length_rounds_as_encode(run):
    """Check the pooled fragment length against ENCODE's rounded mean, run by the image's python."""
    lengths = [fraglen(run, "YCR", f"rep{i}") for i in (1, 2)]
    out = encode(
        "rounded_mean/YCR",
        [],
        f"python -c 'print(int(round(sum([{lengths[0]},{lengths[1]}])/float(2))))' > mean.txt",
    )
    assert fraglen(run, "YCR", "pooled") == int(lines(out / "mean.txt")[0])


@pytest.mark.parametrize("sample", ["rep1", "pooled"])
def test_spp_peaks_match_encode(run, sample):
    """Check SPP regionPeaks against ENCODE's call on the same reads and control."""
    ta = run / f"YCR/tagalign/YCR.{sample}.tagAlign.gz"
    control = control_of(run, "YCR", sample)
    out = encode(
        f"spp/YCR.{sample}",
        [ta, control, CHRSZ],
        f"python3 $(which encode_task_spp.py) {ta.name} {control.name} --chrsz {CHRSZ.name} "
        f"--fraglen {fraglen(run, 'YCR', sample)} --cap-num-peak 300000 --fdr-thresh 0.01 --nth 2",
    )
    ours = regions(run / f"YCR/spp/YCR.{sample}.regionPeak.gz")
    theirs = regions(next(out.glob("*.regionPeak.gz")))
    assert jaccard(ours, theirs) >= 0.99


@pytest.mark.parametrize(("exp", "control"), [("YCR_h", True), ("BGT1", False)])
def test_macs3_peaks_match_encode(run, exp, control):
    """Check MACS3 peak regions against ENCODE's MACS2 call on the same reads."""
    ta = run / f"{exp}/tagalign/{exp}.rep1.tagAlign.gz"
    ctl = [control_of(run, exp, "rep1")] if control else []
    out = encode(
        f"macs2/{exp}.rep1",
        [ta, *ctl, CHRSZ],
        f"python3 $(which encode_task_macs2_chip.py) {ta.name} {' '.join(c.name for c in ctl)} "
        f"--gensz hs --chrsz {CHRSZ.name} --fraglen {fraglen(run, exp, 'rep1')} "
        "--cap-num-peak 500000 --pval-thresh 0.01 --mem-gb 4",
    )
    ours = regions(run / f"{exp}/macs3/{exp}.rep1.narrowPeak.gz")
    theirs = regions(next(out.glob("*.narrowPeak.gz")))
    assert jaccard(ours, theirs) >= 0.999


@pytest.mark.parametrize(("exp", "control"), [("YCR", True), ("BGT1", False)])
def test_signal_tracks_match_encode(run, exp, control):
    """Check fold change and p-value bigWigs against ENCODE's MACS2 tracks, within 2%."""
    ta = run / f"{exp}/tagalign/{exp}.rep1.tagAlign.gz"
    ctl = [control_of(run, exp, "rep1")] if control else []
    out = encode(
        f"signal/{exp}.rep1",
        [ta, *ctl, CHRSZ],
        f"python3 $(which encode_task_macs2_signal_track_chip.py) {ta.name} "
        f"{' '.join(c.name for c in ctl)} --gensz hs --chrsz {CHRSZ.name} "
        f"--fraglen {fraglen(run, exp, 'rep1')} --pval-thresh 0.01 --mem-gb 4",
    )
    for track in ("fc", "pval"):
        theirs = bigwig_stats(next(out.glob(f"*.{track}.signal.bigwig")))
        ours = bigwig_stats(run / f"{exp}/signal/{exp}.rep1.{track}.bigWig")
        assert ours["mean"] == pytest.approx(theirs["mean"], rel=0.02)
        assert ours["max"] == pytest.approx(theirs["max"], rel=0.02)


@pytest.mark.parametrize(
    ("exp", "sample", "caller", "ptype"),
    [
        ("YCR", "rep1", "spp", "regionPeak"),
        ("YCR", "pooled-pr1", "spp", "regionPeak"),
        ("YCR_h", "rep1", "macs3", "narrowPeak"),
        ("BGT1", "rep1", "macs3", "narrowPeak"),
    ],
)
def test_filtered_peaks_and_frip_match_encode(run, exp, sample, caller, ptype):
    """Check filtered peaks, shifted FRiP, peak count and widths against ENCODE's post_call_peak."""
    peaks = run / f"{exp}/{caller}/{exp}.{sample}.{ptype}.gz"
    ta = run / f"{exp}/tagalign/{exp}.{sample}.tagAlign.gz"
    out = encode(
        f"post_call_peak/{exp}.{sample}",
        [peaks, ta, CHRSZ, BLACKLIST],
        f"python3 $(which encode_task_post_call_peak_chip.py) {peaks.name} --ta {ta.name} "
        f"--peak-type {ptype} --fraglen {fraglen(run, exp, sample.split('-')[0])} "
        f"--chrsz {CHRSZ.name} --blacklist {BLACKLIST.name} {REGEX} --mem-gb 4",
    )
    ours = lines(run / f"{exp}/{caller}/{exp}.{sample}.bfilt.{ptype}.gz")
    assert first_difference(ours, lines(next(out.glob(f"*.bfilt.{ptype}.gz")))) is None
    theirs = float(lines(next(out.glob("*.frip.qc")))[0])
    assert frip(run, exp, sample, caller) == pytest.approx(theirs, rel=1e-12)
    assert len(ours) == int(lines(next(out.glob("*.num_peak.qc")))[0])
    if sample in ("rep1", "pooled"):
        sizes = [
            float(x.split("\t")[1])
            for x in lines(next(out.glob("*.peak_region_size.qc")))
        ]
        qc = pl.read_csv(run / f"{exp}/peak_qc/{exp}.tsv", separator="\t")
        row = qc.row(by_predicate=pl.col("label") == sample, named=True)
        ours = [
            row[k]
            for k in ("min_size", "25_pct", "50_pct", "75_pct", "max_size", "mean")
        ]
        assert ours == pytest.approx(sizes, rel=1e-9)


PAIRS = [
    ("YCR", "rep1_vs_rep2"),
    ("YCR", "rep1-pr1_vs_rep1-pr2"),
    ("YCR", "pooled-pr1_vs_pooled-pr2"),
    ("BGT", "rep1_vs_rep2"),
    ("YCR_h", "rep1_vs_rep2"),
    ("BGT1", "rep1-pr1_vs_rep1-pr2"),
]


def pair_files(run, exp, pair):
    """Return a pair's raw peaks, oracle peaks, scoring tagAlign and fragment length.

    Parameters
    ----------
    run : pathlib.Path
        Workflow run directory.
    exp : str
        ChIP experiment.
    pair : str
        Pair label.

    Returns
    -------
    tuple
        Peaks 1, peaks 2, oracle peaks, oracle tagAlign, fragment length.
    """
    caller, ptype = (
        ("spp", "regionPeak") if exp in ("YCR", "BGT") else ("macs3", "narrowPeak")
    )
    first, second = pair.split("_vs_")
    oracle = first.split("-")[0] if "-" in first else "pooled"
    peak = lambda s: run / f"{exp}/{caller}/{exp}.{s}.{ptype}.gz"
    ta = run / f"{exp}/tagalign/{exp}.{oracle}.tagAlign.gz"
    return (
        peak(first),
        peak(second),
        peak(oracle),
        ta,
        fraglen(run, exp, oracle),
        oracle,
    )


@pytest.mark.parametrize("method", ["overlap", "idr"])
@pytest.mark.parametrize(("exp", "pair"), PAIRS)
def test_pair_peaks_and_frip_match_encode(run, exp, pair, method):
    """Check overlap and IDR peaks, filtered, and their FRiP against ENCODE's."""
    if method == "idr" and exp not in ("YCR", "BGT"):
        pytest.skip("IDR runs on tf experiments only")
    p1, p2, oracle_peaks, ta, length, oracle = pair_files(run, exp, pair)
    ptype = p1.name.split(".")[-2]
    flags = (
        "--idr-thresh 0.05 --idr-rank signal.value"
        if method == "idr"
        else "--nonamecheck"
    )
    out = encode(
        f"{method}/{exp}.{pair}",
        [p1, p2, oracle_peaks, ta, CHRSZ, BLACKLIST],
        f"python3 $(which encode_task_{method}.py) {p1.name} {p2.name} {oracle_peaks.name} "
        f"--prefix P --peak-type {ptype} {flags} --fraglen {length} --chrsz {CHRSZ.name} "
        f"--blacklist {BLACKLIST.name} {REGEX} --ta {ta.name}",
    )
    ours = lines(run / f"{exp}/{method}/{exp}.{pair}.bfilt.{ptype}.gz")
    theirs = lines(next(out.glob(f"P.*bfilt.{ptype}.gz")))
    assert first_difference(ours, theirs) is None
    theirs = float(lines(next(out.glob("P.*bfilt.frip.qc")))[0])
    assert frip(run, exp, oracle, f"{method}.{pair}") == pytest.approx(
        theirs, rel=1e-12
    )


@pytest.mark.parametrize(
    ("exp", "method"),
    [("YCR", "idr"), ("YCR", "overlap"), ("YCR_h", "overlap"), ("BGT1", "overlap")],
)
def test_reproducibility_matches_encode(run, exp, method):
    """Check the reproducibility QC and the optimal and conservative sets against ENCODE's."""
    ptype = "regionPeak" if exp == "YCR" else "narrowPeak"
    files = lambda pairs: [
        run / f"{exp}/{method}/{exp}.{p}.bfilt.{ptype}.gz" for p in pairs
    ]
    reps = 1 if exp == "BGT1" else 2
    true = files(["rep1_vs_rep2"] if reps == 2 else [])
    self_pseudo = files([f"rep{i}-pr1_vs_rep{i}-pr2" for i in range(1, reps + 1)])
    pooled = files(["pooled-pr1_vs_pooled-pr2"] if reps == 2 else [])
    ppr = f"--peak-ppr {pooled[0].name}" if pooled else ""
    out = encode(
        f"reproducibility/{exp}.{method}",
        true + self_pseudo + pooled,
        f"python3 $(which encode_task_reproducibility.py) {' '.join(p.name for p in true)} "
        f"--peaks-pr {' '.join(p.name for p in self_pseudo)} {ppr} --prefix M --peak-type {ptype}",
    )
    header, values = lines(out / "M.reproducibility.qc")
    theirs = dict(zip(header.split("\t"), values.split("\t")))
    ours = pl.read_csv(
        run / f"{exp}/reproducibility/{exp}.{method}.{ptype}.tsv", separator="\t"
    ).row(0, named=True)
    assert list(ours) == list(theirs)
    for key, value in theirs.items():
        if key.endswith("ratio"):
            assert ours[key] == pytest.approx(float(value), rel=1e-12)
        else:
            assert str(ours[key]) == value
    for kind in ("optimal", "conservative"):
        ours = lines(run / f"{exp}/reproducibility/{exp}.{method}.{kind}.{ptype}.gz")
        assert first_difference(ours, lines(out / f"M.{kind}_peak.{ptype}.gz")) is None


@pytest.mark.parametrize(("exp", "control"), [("YCR", "YCR_ctl"), ("BGT1", None)])
def test_jsd_matches_encode(run, exp, control):
    """Check the fingerprint metrics of each replicate against ENCODE's, within 1%."""
    bams = sorted((run / f"{exp}/filter").glob(f"{exp}.rep[0-9].bam"))
    ctl = [run / f"{control}/filter/{control}.ctl1.bam"] if control else []
    out = encode(
        f"jsd/{exp}",
        [*bams, *ctl, BLACKLIST],
        f"python3 $(which encode_task_jsd.py) {' '.join(b.name for b in bams)} "
        f"{'--ctl-bam ' + ctl[0].name if ctl else ''} --blacklist {BLACKLIST.name} "
        "--mapq-thresh 30 --nth 4",
    )
    ours = {
        row[0]: [float(x) for x in row[1:]]
        for row in (x.split("\t") for x in lines(run / f"{exp}/jsd/{exp}.tsv")[1:])
    }
    for i in range(1, len(bams) + 1):
        theirs = [
            float(x) for x in lines(next(out.glob(f"rep{i}.*.jsd.qc")))[0].split("\t")
        ]
        assert ours[f"rep{i}"] == pytest.approx(theirs, rel=0.01)


def test_gc_bias_matches_encode(run):
    """Check Picard's normalized coverage per GC bin against ENCODE's run."""
    bam = run / "YCR/filter/YCR.rep1.bam"
    out = encode(
        "gc_bias/YCR.rep1",
        [bam, GENOME],
        f"python3 $(which encode_task_gc_bias.py) --nodup-bam {bam.name} "
        f"--ref-fa {GENOME.name} --picard-java-heap 4G",
    )
    column = "NORMALIZED_COVERAGE"

    def coverage(path):
        rows = [x.split("\t") for x in lines(path) if x and not x.startswith("#")]
        i = rows[0].index(column)
        return [float(r[i]) for r in rows[1:]]

    ours = coverage(run / "YCR/gc_bias/YCR.rep1.gc.txt")
    assert ours == pytest.approx(coverage(next(out.glob("*.gc.txt"))), rel=1e-6)


@pytest.mark.parametrize("exp", ["YCR", "BGT1"])
def test_count_signal_tracks_match_encode(run, exp):
    """Check the 5' count bigWigs against ENCODE's count_signal_track."""
    ta = run / f"{exp}/tagalign/{exp}.rep1.tagAlign.gz"
    out = encode(
        f"count_signal/{exp}.rep1",
        [ta, CHRSZ],
        f"python3 $(which encode_task_count_signal_track.py) {ta.name} --chrsz {CHRSZ.name}",
    )
    for strand in ("positive", "negative"):
        theirs = bigwig_stats(next(out.glob(f"*.{strand}.bigwig")))
        assert bigwig_stats(run / f"{exp}/signal/{exp}.rep1.{strand}.bigWig") == theirs
