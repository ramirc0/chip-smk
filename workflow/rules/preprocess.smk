# Shell-rule stubs. Each shows the directive set: software, tee log, benchmark, :q.


# Run-independent resource: built once under results/refs/, reused across runs.
rule reference:
    input:
        REFERENCE,
    output:
        f"{REFDIR}/reference.index",
    log:
        f"{LOGDIR}/reference.txt",
    benchmark:
        f"{BENCHDIR}/reference.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        wc -l {input:q} > {output:q}
        """


# Per-sample input staged to the canonical path. Swap cp for real preprocessing.
rule prepare:
    input:
        lambda w: INPUT_OF[w.sample],
    output:
        f"{OUTDIR}/{{sample}}/{{sample}}.input.tsv",
    log:
        f"{LOGDIR}/prepare/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/prepare/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        cp {input:q} {output:q}
        """
