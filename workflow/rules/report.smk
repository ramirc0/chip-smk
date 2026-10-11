# Aggregation stub. Shows: expand() fan-in over SAMPLES.


rule summary:
    input:
        metrics=expand(
            f"{OUTDIR}/{{sample}}/{{sample}}.metrics.tsv", sample=SAMPLES
        ),
    output:
        plot=f"{OUTDIR}/report/summary.svg",
        plot_png=f"{OUTDIR}/report/summary.png",
    log:
        f"{LOGDIR}/summary.txt",
    benchmark:
        f"{BENCHDIR}/summary.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/plot_summary.py \
            --input {input.metrics:q} \
            --output {output.plot:q}
        """
