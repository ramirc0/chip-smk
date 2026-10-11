# CLAUDE.md

Skeleton Snakemake workflow. Rules are stubs that demonstrate the directive set.
README.md holds the layout, idioms, and style conventions. Follow it.

## Commands

```bash
pixi run lint                                # ruff check, ruff format --check, numpydoc lint
pixi run test                                # full run in .tests/integration
pixi run snakemake -n -p --profile profiles/local
pixi run snakemake --profile profiles/local --cores 2
```

## Two pixi workspaces

- Root `pixi.toml` is the launcher: Snakemake and its plugins, plus the `lint` env.
- `workflow/envs/pixi.toml` is the rule env. Rule dependencies go there, never in the root.
- `pyproject.toml` is ruff and numpydoc config only. pixi ignores it.
- Rules use `locked=True`. After a manifest edit, run `pixi lock` in `workflow/envs/`.
- Any manifest or lock change reruns every rule.
- `SOFTWARE_ENV` resolves `../envs` relative to the rule file. New rule files MUST sit in `workflow/rules/`.

## Gotchas

- Shell blocks MUST call scripts as `{workflow.basedir:q}/scripts/x.py`. A cwd-relative path breaks `--directory`.
- `.smk` files are not linted. Keep helper docstrings numpydoc by hand.
- Plots use poikilos `plain`. `pk.save_figure()` writes SVG and PNG. Rules MUST declare both outputs.
- `MPLBACKEND=Agg` is set in the rule env activation. Scripts rely on it.
- Resources MUST be flat values. Never derive `mem_mb` or `runtime` from input size.
- The `* attempt` retry scaling in `profiles/slurm/config.yaml` is intentional.
- Auto mode may block `pixi lock`, `pixi install`, or workflow runs that fetch the poikilos git dependency. Ask the user to run them.

## Verifying a change

`pixi run test` runs every rule on the toy data in `.tests/integration/`.
New rules or inputs MUST get matching test inputs there.
