# CLAUDE.md

ChIP-seq Snakemake workflow on the snakemake-template layout.
README.md holds the layout, steps, outputs, and style conventions. Follow it.

## Commands

```bash
pixi run lint                                # ruff check, ruff format --check, numpydoc lint
pixi run snakemake -n -p --profile profiles/local
pixi run snakemake --profile profiles/local
TMPDIR=/data/... pixi run test               # chr19 fixture vs ENCODE's task scripts, pytest-xdist
pixi run snakemake --profile profiles/slurm --configfile .tests/validation/config.yaml  # ENCODE validation run
```

## Two pixi workspaces

- Root `pixi.toml` is the launcher: Snakemake and its plugins, plus the `lint` env.
- `workflow/envs/pixi.toml` holds the rule envs `default`, `idr`, `spp` and `test`. Rule dependencies go there, never in the root.
- `pyproject.toml` is ruff and numpydoc config only. pixi ignores it.
- Rules use `locked=True`. After a manifest edit, run `pixi lock` in `workflow/envs/`.
- Any manifest or lock change reruns every rule.
- `SOFTWARE_ENV`, `IDR_ENV` and `SPP_ENV` resolve `../envs` relative to the rule file. New rule files MUST sit in `workflow/rules/`.

## Match ENCODE

- Every step reproduces ENCODE ChIP-seq pipeline v2.2.2 (`~/Git_Repos/chip-seq-pipeline2`). README.md "Differences from ENCODE" lists the deliberate gaps.
- `.tests/integration/test_encode.py` checks each step against ENCODE's own task script on the same inputs. A step change MUST keep its test passing.
- A bug fix starts with a test that fails on the old code.
- Correctness beats parity. Where ENCODE has a known bug, fix it and list the gap in README.md.
- ENCODE's task scripts crash without the flags its WDL always passes, such as an empty `--filter-chrs`.

## Gotchas

- Names never contain dots. File names are `<experiment>.<unit>...` and the wildcards rely on it.
- Shuffles (pseudoreplicates, subsamples, xcor) seed OpenSSL with the tagAlign's byte size. Keep `-md md5`, keep `shuf` reading a pipe, and keep tagAlign text byte-identical.
- `paste - -` joins mates. It equals ENCODE's `sed 'N;s/\n/\t/'` byte for byte.
- Every tagAlign has a `.counts` file (`lines<TAB>bytes`), written by the awk that writes it. Never size a tagAlign by reading it again.
- bowtie2 runs with `--reorder`; equal-coordinate order and every shuffle depend on it.
- `choose_ctl` is a checkpoint. Dry runs stop at it until the tagAligns exist.
- bedtools drops `-names` with a single `-b` file, so FRiP tags each peak with its set instead.
- cutadapt cannot sniff the compression of a pipe; feed it decompressed lanes.
- `rapidgzip` has no Python 3.9 build, so the `idr` env uses `pigz`.
- chrM stays in BAMs and tagAligns, as ENCODE's `filter_chrs` is empty.
- Plots use poikilos `plain`. `pk.save_figure()` writes SVG and PNG. Rules MUST declare both outputs.
- Resources MUST be flat values. Never derive `mem_mb` or `runtime` from input size.
- run_spp.R's snow workers detach, so benchmarks and sacct miss them. Size `spp` from the job cgroup's `memory.peak`.
- The signal sort and de-overlap awk are needed: slop clamps intervals past a contig end into ties.
- The `* attempt` retry scaling in `profiles/slurm/config.yaml` is intentional.
- CLI `--set-threads` replaces the profile's whole mapping. Pass every rule.

## Verifying a change

1. `pixi run lint` and `pixi run test`. The login node's `/` is nearly full: set `TMPDIR` to `/data` scratch.
2. For a big change, run the validation on SLURM, then `.tests/validation/compare_encode.py`.
