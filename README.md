# chip-smk

A [Snakemake](https://github.com/snakemake/snakemake) pipeline for ChIP-seq. It
reproduces the [ENCODE ChIP-seq pipeline](https://github.com/ENCODE-DCC/chip-seq-pipeline2)
v2.2.2 step by step, for transcription factors (SPP peaks, IDR and naive overlap)
and histone marks (MACS3 peaks, naive overlap). Replicates, controls, pools and
pseudoreplicates are wired as ENCODE wires them. The work ENCODE repeats or throws
away is gone.

## Layout

```
config/
  config.yaml.template   parameters + paths; copy to config.yaml (gitignored)
  samples.example.tsv    one row per FASTQ (pair): experiment, replicate, type, ...
workflow/
  Snakefile              includes the rule files, defines rule all
  rules/
    common.smk           config, sample sheet, naming helpers
    preprocess.smk       indexes, fastp, crop + alignment, filter, tagAlign, subsample, pool, pseudoreplicates
    xcor.smk             cross-correlation reads, run_spp.R, fragment lengths
    peaks.smk            control choice, SPP, MACS3, exclusion filter, overlap, IDR, peak sets
    signal.smk           fold-change, p-value and count bigWigs
    qc.smk               PBC/NRF, FRiP, JSD, GC bias, peak widths
    report.smk           run-level QC table and figure
  scripts/               argparse scripts, each with build_parser() + main()
  envs/pixi.toml         rule envs (pixi workspace): default, idr, spp and test
.tests/
  integration/           chr19 fixture run, step by step against ENCODE's task scripts
  validation/            four ENCODE experiments on GRCh38, compared with ENCODE's QC
profiles/
  local/config.yaml      local execution
  slurm/config.yaml      SLURM (per-rule resources, retried with more memory)
pixi.toml                launcher env, lint and test tasks
pyproject.toml           ruff + numpydoc config
```

## Running

```bash
cp config/config.yaml.template config/config.yaml  # then edit
cp config/samples.example.tsv config/samples.tsv   # then edit
pixi shell                                         # launcher env
snakemake -n -p --profile profiles/local           # dry run
snakemake --profile profiles/local                 # local
snakemake --profile profiles/slurm                 # SLURM
```

Snakemake MUST run from the repository root. Control choice is a checkpoint:
which control file each peak call reads depends on read depths. A dry run lists
the jobs downstream of it only once the tagAligns exist.

The SLURM profile holds Weng Lab cluster settings, such as partitions and the
`cascadelake` constraint. Adapt them for another cluster with the
[SLURM executor docs](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).
Tune `set-resources` from the `benchmarks/<run_id>/` TSVs after a first run.
`run_spp.R` (rules `spp` and `xcor`) starts snow workers that detach from the
job, so the benchmarks and the SLURM efficiency report miss their CPU and memory.
Size `spp` from the job cgroup's `memory.peak` instead.

## Testing

Both tests use data on the Weng Lab cluster, so they run only there.

`pixi run test` runs the integration tests. The fixture holds the chr19 and chrM
reads of two ENCODE experiments and their controls: ENCSR014YCR (ATF7,
paired-end) and ENCSR000BGT (BATF, single-end, three controls). The fixture sheet
reuses them as histone experiments too, one with a shared control and one with
none. The tests run the workflow on it. Then they run ENCODE's own task scripts,
from the `encodedcc/chip-seq-pipeline:v2.2.2` image, on the same inputs, and
compare. Build the fixture once and pull the image once:

```bash
pixi run --manifest-path workflow/envs/pixi.toml -e test bash .tests/integration/make_fixture.sh
apptainer pull /data/zusers/ramirezc/chip-smk-refs/containers/chip-seq-pipeline_v2.2.2.sif \
    docker://encodedcc/chip-seq-pipeline:v2.2.2   # or set ENCODE_SIF
TMPDIR=/data/... pixi run test
```

The validation runs four ENCODE4 experiments on full GRCh38, then prints their
QC beside ENCODE's published values:

| Experiment | Target | Reads | Control |
|---|---|---|---|
| ENCSR014YCR | ATF7 (tf), GM12878 | PE 100 | ENCSR398JTO, 2 replicates |
| ENCSR000BGT | BATF (tf), GM12878 | SE 36 | ENCSR000BGH, 3 replicates of uneven depth |
| ENCSR736ZEG | H3K27ac (histone), keratinocyte | PE 101 | ENCSR982AYN, 2 replicates |
| ENCSR243BPI | ZNF219 (tf), HepG2 | SE 50 | ENCSR346XJA, 1 replicate |

```bash
pixi run --manifest-path workflow/envs/pixi.toml -e test bash .tests/validation/make_inputs.sh
snakemake --profile profiles/slurm --configfile .tests/validation/config.yaml
pixi run --manifest-path workflow/envs/pixi.toml python .tests/validation/compare_encode.py
```

## Inputs

The sample sheet is a TSV with one row per FASTQ or FASTQ pair. Lines starting
with `#` are comments.

| Column | Required | Meaning |
|---|---|---|
| `experiment` | yes | Experiment name: letters, digits, `_` and `-`. No dots. |
| `replicate` | yes | 1 to N within the experiment. Rows sharing both are lanes of one library. |
| `type` | yes | `tf`, `histone` or `control`. |
| `control` | tf, optional for histone | The `control` experiment to call peaks against. Several experiments MAY share one. |
| `crop_length` | no | Crop reads to this length before alignment, as ENCODE4 did. Empty means no crop. |
| `r1` | yes | Gzipped FASTQ. |
| `r2` | no | Gzipped mate FASTQ. Empty means single-end. |

`type`, `control` and `crop_length` MUST agree on every row of an experiment.
A control experiment is aligned once, however many experiments use it. Its
replicates become ENCODE's `ctl1` to `ctlN`.

References: the GRCh38 no-alt analysis set FASTA
(`GRCh38_no_alt_analysis_set_GCA_000001405.15`) and ENCODE's exclusion list
ENCFF356LFX.

## Config

`config/config.yaml.template` lists every key with ENCODE's defaults. Paths:
`samples`, `outdir`, `run_id`, `tmp_dir` and `references`. Parameters, under
`params`: the aligner (`bowtie2` or `bwa`), MAPQ, subsampling, xcor, control
choice, peak thresholds and caps, `true_rep_only` and `count_signal_track`.

## Steps

Preprocessing, in `preprocess.smk`:

- `fastp_qc` writes fastp's report per lane. It trims nothing.
- `bowtie2_build` and `bwa_index` index the genome once, under `results/refs/`.
- `align` crops when `crop_length` is set, aligns with bowtie2 (`--mm`, `-X 2000`
  for pairs) or bwa (`aln`, or `mem` for long pairs), and sorts.
- `filter` keeps MAPQ ≥ 30 primary reads, proper pairs after `fixmate`, marks
  duplicates with sambamba and removes them.
- `tagalign` writes ENCODE's tagAlign and a `lines<TAB>bytes` counts file.
- `subsample_ta` shuffles with ENCODE's seeded keystream and keeps N reads.
- `pool_ta` concatenates replicates. `pseudoreps` splits a replicate in two.

Fragment length, in `xcor.smk`:

- `xcor_tagalign` aligns R1, trimmed to 50 bp unless a read is shorter, and keeps
  MAPQ ≥ 30 reads without duplicate removal.
- `xcor` subsamples 15M non-chrM reads and runs `run_spp.R`. `plot_xcor` draws the
  profile. `fraglen_mean` gives pools the replicates' rounded mean.

Peaks, in `peaks.smk`:

- `choose_ctl` picks each replicate's control (one, its own, or the pool) and its
  subsample depth, as ENCODE does.
- `spp` calls tf peaks. `macs3` calls histone peaks and writes the bedGraphs for
  signal. `macs3_pseudorep` calls pseudoreplicate peaks.
- `blacklist_filter` drops peaks that touch the exclusion list or lie off chr1-22,
  X and Y.
- `overlap` and `idr` compare true replicates, each replicate's pseudoreplicates,
  and the pooled pseudoreplicates. `reproducibility` picks the optimal and
  conservative sets and writes ENCODE's reproducibility QC.

Signal, in `signal.smk`: `signal` writes fold-change and p-value bigWigs per
replicate and pool. `count_signal` writes 5' read counts per strand.

QC, in `qc.smk` and `report.smk`: `pbc_qc` (NRF, PBC1, PBC2), `frip` (shifted
FRiP of every peak set a tagAlign scores, in one pass), `jsd` (deeptools
fingerprints), `gc_bias` (Picard) and `peak_qc` (peak counts and widths).
`qc_summary` gathers everything into `qc.tsv`, keyed like ENCODE's `qc.json`.

## Outputs

Each experiment writes to `results/<run_id>/<experiment>/<step>/`, with files
named `<experiment>.<unit>...`. Units follow ENCODE's labels: `rep1`, `ctl1`,
`pooled`, `rep1-pr1`, `rep1_vs_rep2`, `pooled-pr1_vs_pooled-pr2`. The indexes and
chromosome sizes go to `results/refs/`. Logs and benchmarks mirror the run under
`logs/<run_id>/` and `benchmarks/<run_id>/`.

| Step directory | Key outputs |
|---|---|
| `fastp/` | fastp HTML and JSON per lane |
| `align/` | Sorted BAM, flagstat JSON |
| `filter/` | Filtered, deduplicated BAM and index; flagstat JSON before and after duplicate removal |
| `tagalign/` | tagAligns of replicates, controls, pools, pseudoreplicates and subsamples |
| `xcor/` | `cc.qc`, fragment length, profile and plot |
| `choose_ctl/` | Control and subsample depth per replicate and pool |
| `spp/`, `macs3/` | Peaks and filtered (`.bfilt.`) peaks per sample |
| `overlap/`, `idr/` | Filtered peaks per pair; IDR scores and plots |
| `reproducibility/` | Optimal and conservative peaks, reproducibility QC |
| `signal/` | Fold-change and p-value bigWigs, optional count bigWigs |
| `frip/`, `pbc_qc/`, `jsd/`, `gc_bias/`, `peak_qc/` | QC tables and plots |

`results/<run_id>/qc.tsv` holds every metric of the run, one per row:
`experiment, section, unit, metric, value`. `qc.svg` plots the key ones.

## Environment

The launcher env (`pixi.toml` at the root) has Snakemake main and its plugins
for pixi deployment and SLURM. The rule envs live in a second pixi workspace,
`workflow/envs/`. Its `pixi.lock` pins every package. `default` runs most rules.
`idr` runs IDR on the Python 3.9 / NumPy 1.19 stack it needs. `spp` runs
`run_spp.R`. `test` adds pytest, pytest-xdist and seqkit. pixi MUST be on `PATH`.
Jobs install the envs on first use.

Add a dependency with `pixi add <pkg>` in `workflow/envs/`. After a hand edit,
run `pixi lock` there (rules use `locked=True`). Any change to the manifest or
lock reruns every rule. So does moving the workdir, because the env hash includes
the workspace's absolute path.

## Conventions

The workflow follows the
[Nextstrain Snakemake style guide](https://docs.nextstrain.org/en/latest/reference/snakemake-style-guide.html)
and [Snakemake's best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html):

- No `run:` blocks. Custom logic lives in `workflow/scripts/`, called from `shell:`.
- Raw, triple-quoted `shell` blocks (`r"""`), one option per line.
- Every interpolation uses `:q` quoting.
- Config reaches `shell` through `params:` (a `lambda`), never directly.
- `input:` and `output:` are literal path strings. All paths are relative.
- Rules tee logs to file and terminal (`exec &> >(tee {log:q})`).
- Every rule has `log:` and `benchmark:`. Each names its pixi env in `software:`.
- Plots use [poikilos](https://github.com/ramirc0/poikilos) with the `plain`
  theme. `pk.save_figure()` writes the SVG and PNG pair. Rules declare both.
- Docstrings follow numpydoc. `pixi run lint` checks the scripts with ruff and
  `numpydoc lint`. Neither tool parses `.smk` files.

## Differences from ENCODE

Outputs match ENCODE's task scripts on the fixture, except where listed here.

- MACS3 3.0.5 replaces MACS2 2.2.4. Peak regions match.
- sambamba marks duplicates, with ENCODE's settings. ENCODE's ChIP default is
  Picard. Duplicate counts come from `samtools flagstat`.
- `samtools flagstat` replaces SAMstats. They agree without secondary alignments.
- cutadapt crops (`--length`, `-m`) instead of Trimmomatic. The reads kept are the same.
- One MACS3 call per replicate and pool gives histone peaks and both signal
  tracks. ENCODE runs it twice.
- The xcor reads stream from FASTQ to tagAlign. ENCODE also writes their BAM,
  marks duplicates and computes QC, then discards all three.
- A control experiment is aligned once. A control subsample is made once per
  depth, shared by every replicate and experiment that uses it. ENCODE repeats
  both per experiment and per replicate.
- bowtie2 runs with `--reorder`, so reads with equal coordinates, and with them
  pseudoreplicates and subsamples, do not depend on thread timing.
- bwa's bad-CIGAR filter matches whole read names. ENCODE's `grep -F` also drops
  reads whose names contain a flagged name.
- IDR on fewer than 20 merged peaks writes empty results, and ratios over zero
  peaks are infinite. ENCODE fails the run.
- A paired-end control subsample of odd size is rounded down to whole pairs.
  ENCODE raises.
- phantompeakqualtools 1.2.2 and spp 1.16 replace 1.2.1 and 1.15.
- bowtie2 2.5.5 replaces 2.3.4.3. The alignments match on the fixture.
- ENCODE's `qc.json` writes the two-read count to `positions_with_one_read`, as
  its parser sets that key twice. `qc.tsv` keeps both counts under their own names.
- `pct_duplicate_reads` is unrounded. ENCODE rounds the duplicate fraction to two
  decimals, then multiplies by 100.
- Not ported: bigBed, starch and hammock files; `qc.html` and `qc.json`; BAM,
  tagAlign and peak entry points; `filter_chrs`, redaction, bowtie2 `--local`,
  custom aligners, a fixed fragment length or seed, `no_dup_removal`, a peak
  caller override, a second exclusion list, and genomes other than GRCh38.

## References

Chen, Shifu. 2025. “Fastp 1.0: An Ultra-Fast All-Round Tool for FASTQ Data Quality Control and Preprocessing.” iMeta 4 (5): e70078. https://doi.org/10.1002/imt2.70078.

Danecek, Petr, James K. Bonfield, Jennifer Liddle, et al. 2021. “Twelve Years of SAMtools and BCFtools.” GigaScience 10 (2): giab008. https://doi.org/10.1093/gigascience/giab008.

Kharchenko, Peter V., Michael Y. Tolstorukov, and Peter J. Park. 2008. “Design and Analysis of ChIP-Seq Experiments for DNA-Binding Proteins.” Nature Biotechnology 26 (12): 1351–59. https://doi.org/10.1038/nbt.1508.

Landt, Stephen G., Georgi K. Marinov, Anshul Kundaje, et al. 2012. “ChIP-Seq Guidelines and Practices of the ENCODE and modENCODE Consortia.” Genome Research 22 (9): 1813–31. https://doi.org/10.1101/gr.136184.111.

Langmead, Ben, and Steven L. Salzberg. 2012. “Fast Gapped-Read Alignment with Bowtie 2.” Nature Methods 9 (4): 357–59. https://doi.org/10.1038/nmeth.1923.

Li, Heng, and Richard Durbin. 2009. “Fast and Accurate Short Read Alignment with Burrows-Wheeler Transform.” Bioinformatics 25 (14): 1754–60. https://doi.org/10.1093/bioinformatics/btp324.

Li, Qunhua, James B. Brown, Haiyan Huang, and Peter J. Bickel. 2011. “Measuring Reproducibility of High-Throughput Experiments.” The Annals of Applied Statistics 5 (3): 1752–79. https://doi.org/10.1214/11-AOAS466.

Martin, Marcel. 2011. “Cutadapt Removes Adapter Sequences from High-Throughput Sequencing Reads.” EMBnet.Journal 17 (1): 10–12. https://doi.org/10.14806/ej.17.1.200.

Mölder, Felix, Kim Philipp Jablonski, Brice Letcher, et al. 2021. “Sustainable Data Analysis with Snakemake.” F1000Research 10 (January): 33. https://doi.org/10.12688/f1000research.29032.3.

Ramírez, Fidel, Devon P. Ryan, Björn Grüning, et al. 2016. “deepTools2: A next Generation Web Server for Deep-Sequencing Data Analysis.” Nucleic Acids Research 44 (W1): W160–65. https://doi.org/10.1093/nar/gkw257.

Tarasov, Artem, Albert J. Vilella, Edwin Cuppen, Isaac J. Nijman, and Pjotr Prins. 2015. “Sambamba: Fast Processing of NGS Alignment Formats.” Bioinformatics 31 (12): 2032–34. https://doi.org/10.1093/bioinformatics/btv098.

Zhang, Yong, Tao Liu, Clifford A. Meyer, et al. 2008. “Model-Based Analysis of ChIP-Seq (MACS).” Genome Biology 9 (9): R137. https://doi.org/10.1186/gb-2008-9-9-r137.
