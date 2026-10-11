# Run-level QC table and figure.


def qc_files():
    """Return every per-step QC file of the run, by kind.

    Returns
    -------
    dict of str to list of str
        Paths for each `qc_summary.py` flag. File names carry the experiment
        and the replicate, pool or pair they describe.
    """
    files = {k: [] for k in ("flagstat", "pbc", "fastp", "xcor", "jsd", "frip", "reproducibility", "peak_qc")}
    for (exp, lib), entry in LIBRARIES.items():
        prefix = f"{OUTDIR}/{exp}"
        files["flagstat"] += [
            f"{prefix}/align/{exp}.{lib}.flagstat.json",
            f"{prefix}/filter/{exp}.{lib}.flagstat.json",
            f"{prefix}/filter/{exp}.{lib}.dup.flagstat.json",
        ]
        files["pbc"].append(f"{prefix}/pbc_qc/{exp}.{lib}.json")
        files["fastp"] += [f"{prefix}/fastp/{exp}.{lib}.{i}.json" for i in range(1, len(entry["r1"]) + 1)]
    for exp in CHIPS:
        prefix = f"{OUTDIR}/{exp}"
        files["xcor"] += [f"{prefix}/xcor/{exp}.{r}.cc.qc" for r in EXPERIMENTS[exp]["libs"]]
        files["jsd"].append(f"{prefix}/jsd/{exp}.tsv")
        files["frip"] += [f"{prefix}/frip/{exp}.{s}.tsv" for s in samples(exp)]
        files["peak_qc"].append(f"{prefix}/peak_qc/{exp}.tsv")
        if not config["params"]["true_rep_only"]:
            files["reproducibility"] += [
                f"{prefix}/reproducibility/{exp}.{m}.{peak_type(exp)}.tsv" for m in methods(exp)
            ]
    return files


rule qc_summary:
    input:
        **qc_files(),
    output:
        table=f"{OUTDIR}/qc.tsv",
        plot=f"{OUTDIR}/qc.svg",
        plot_png=f"{OUTDIR}/qc.png",
    log:
        f"{LOGDIR}/qc_summary.txt",
    benchmark:
        f"{BENCHDIR}/qc_summary.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/qc_summary.py \
            --flagstat {input.flagstat:q} \
            --pbc {input.pbc:q} \
            --fastp {input.fastp:q} \
            --xcor {input.xcor:q} \
            --jsd {input.jsd:q} \
            --frip {input.frip:q} \
            --reproducibility {input.reproducibility:q} \
            --peak-qc {input.peak_qc:q} \
            --table {output.table:q} \
            --plot {output.plot:q}
        """
