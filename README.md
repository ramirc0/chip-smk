# snakemake-template

A skeleton Snakemake workflow: a blank canvas with the layout and conventions
pre-wired. The rules are stubs (`cp`, `wc -l`, a script call) whose only job is
to demonstrate the directive set; replace their bodies with real work. The
example wires `reference` + `prepare -> analyze -> summary` over a sample sheet.

## Layout

```
.tests/integration/      config + toy inputs for the end-to-end test
config/
  config.yaml.template   parameters + paths; copy to config.yaml (gitignored)
  samples.example.tsv    one row per sample: sample_id, input, [label]
workflow/
  Snakefile              includes the rules + rule all
  rules/
    common.smk           config, sample sheet, path constants, helpers
    preprocess.smk       shell-only prep (shared reference + per-sample)
    analyze.smk          main per-sample compute (flag-driven script)
    report.smk           run-level aggregation
  scripts/               argparse scripts, each with build_parser() + main()
  envs/pixi.toml         rule env (pixi workspace)
profiles/
  local/config.yaml      local execution
  slurm/config.yaml      SLURM (per-rule resources, retried with more memory)
pixi.toml                launcher env + test and lint tasks
pyproject.toml           ruff + numpydoc config
```

## Idioms

- **Outputs grouped by `run_id`** under `results/<run_id>/<sample>/`; logs and
  benchmarks mirror it (`logs/<run_id>/`, `benchmarks/<run_id>/`). Override
  `run_id` per run to keep runs separate. Run-independent artifacts go in
  `results/refs/`.
- **Sample sheet** is a TSV read once in `common.smk`; optional columns are
  looked up with `_column()` (empty cell -> `None`).
- **Config -> flags**: rules stay declarative; helpers like `analyze_flags()`
  turn `config` into a script flag string.
- **Every rule** carries `log:`, `benchmark:`, and `software:`.
- **Scripts** expose `build_parser()` so defaults are introspectable.
- **Plots** use [poikilos](https://github.com/ramirc0/poikilos) with the
  `plain` theme: `pk.use("plain")` before drawing, `pk.despine()` on every axes,
  `pk.save_figure()` for the SVG + PNG pair. Rules declare both files.

## Conventions

Follow the [Nextstrain Snakemake style guide][style-guide] and
[Snakemake's best-practices guide][best-practices]. The stubs already apply it:

- No `run:` blocks; custom logic lives in `workflow/scripts/` called from `shell:`.
- Raw, triple-quoted `shell` blocks (`r"""`), one option per line.
- Every interpolation quoted with `:q`; multi-value params passed as lists.
- Config mapped into `shell` via `params:` (a `lambda`), never interpolated directly.
- `config[key]` for required keys, `config.get(key, default)` / `key in config`
  for optional ones; never bare `config.get(key)`.
- `input:`/`output:` are literal path strings (no `rules.x.output` refs); all
  paths relative.
- Logs tee'd to file and terminal (`exec &> >(tee {log:q})`); run with
  `--show-failed-logs` (set in the profiles).
- Every rule has a `benchmark:`; no `message:` attribute.
- Scripts are called as `{workflow.basedir:q}/scripts/x.py` so the workflow
  runs from any `--directory`.
- Docstrings follow numpydoc. `pixi run lint` checks the scripts with ruff
  (pydocstyle, numpy convention) and `numpydoc lint`. Neither tool parses
  `.smk` files. Keep the rule helpers in the same style by hand.

[style-guide]: https://docs.nextstrain.org/en/latest/reference/snakemake-style-guide.html
[best-practices]: https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html

## Running

```bash
cp config/config.yaml.template config/config.yaml  # then edit
cp config/samples.example.tsv config/samples.tsv   # then edit
pixi shell                                         # launcher env
snakemake --profile profiles/local                 # local
snakemake --profile profiles/slurm                 # SLURM
snakemake -n -p --profile profiles/local           # dry run
pixi run test                                      # end-to-end test
```

Tune `set-resources` in the profiles from the `benchmarks/<run_id>/` TSVs after
a first run.

## Environment

The launcher env (`pixi.toml` at the root) has Snakemake main, the pixi
software-deployment plugin and the SLURM executor plugins. The rule env is a
second pixi workspace, `workflow/envs/`. Its `pixi.lock` pins every package.
pixi must be on `PATH`. Jobs install the env on first use.

Add a dependency with `pixi add <pkg>` in `workflow/envs/`. After a hand edit,
run `pixi lock` there (rules use `locked=True`). Any change to the manifest or
lock reruns every rule. So does moving the workdir, because the env hash
includes the workspace's absolute path.
