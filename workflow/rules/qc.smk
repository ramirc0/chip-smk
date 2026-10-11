# Library, enrichment and peak QC.


# ENCODE's library complexity, on the filtered reads before duplicate removal.
rule pbc_qc:
    input:
        f"{OUTDIR}/{{exp}}/filter/{{exp}}.{{lib}}.filt.bam",
    output:
        f"{OUTDIR}/{{exp}}/pbc_qc/{{exp}}.{{lib}}.json",
    params:
        mito=lambda w: config["mito_chr_name"],
        bedpe=lambda w: ["-bedpe"] if LIBRARIES[w.exp, w.lib]["layout"] == "pe" else [],
        paired=lambda w: ["--paired-end"] if LIBRARIES[w.exp, w.lib]["layout"] == "pe" else [],
    log:
        f"{LOGDIR}/pbc_qc/{{exp}}.{{lib}}.txt",
    benchmark:
        f"{BENCHDIR}/pbc_qc/{{exp}}.{{lib}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        export POLARS_MAX_THREADS={threads}
        bedtools bamtobed {params.bedpe:q} -i {input:q} \
            | python {workflow.basedir:q}/scripts/pbc_qc.py \
                {params.paired:q} \
                --mito {params.mito:q} \
                --output {output:q}
        """


# ENCODE's shifted FRiP: each read moved half a fragment downstream, then the reads
# that touch a peak, for every peak set the tagAlign scores, in one pass. bedtools
# slop truncates ENCODE's fractional half, and its clamps at contig ends are kept.
rule frip:
    input:
        ta=lambda w: ta(w.exp, sample_ta(w.exp, w.sample)),
        counts=lambda w: counts(w.exp, sample_ta(w.exp, w.sample)),
        fraglen=lambda w: fraglen(w.exp, w.sample),
        chromsizes=CHROMSIZES,
        peaks=lambda w: list(frip_sets(w.exp, w.sample).values()),
    output:
        f"{OUTDIR}/{{exp}}/frip/{{exp}}.{{sample}}.tsv",
    params:
        labels=lambda w: list(frip_sets(w.exp, w.sample)),
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/frip/{{exp}}.{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/frip/{{exp}}.{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        read -r reads _ < {input.counts:q}
        half=$(( ($(< {input.fraglen:q}) + 1) / 2 ))
        labels=({params.labels:q})
        peaks=({input.peaks:q})
        # One file of every set's regions, each tagged with its label: bedtools drops
        # -names with a single -b file, and mishandles gzipped -b files.
        for i in "${{!labels[@]}}"; do
            pigz -dc "${{peaks[i]}}" | awk -v l="${{labels[i]}}" 'BEGIN {{ OFS = "\t" }} {{ print $1, $2, $3, l }}'
        done > "$tmp/peaks.bed"
        if [[ -s "$tmp/peaks.bed" ]]; then
            rapidgzip -d -c -P {threads} {input.ta:q} \
                | bedtools slop -i stdin -g {input.chromsizes:q} -s -l "-$half" -r "$half" \
                | awk 'BEGIN {{ OFS = "\t" }} $2 >= 0 && $3 >= 0 && $2 <= $3 {{ print $1, $2, $3, NR }}' \
                | bedtools intersect -nonamecheck -wa -wb -a stdin -b "$tmp/peaks.bed" \
                | awk '$4 != id {{ id = $4; delete seen }} !seen[$8]++ {{ hits[$8]++ }} END {{ for (l in hits) print hits[l], l }}' \
                > "$tmp/hits"
        else
            : > "$tmp/hits"
        fi
        awk -v reads="$reads" 'BEGIN {{ OFS = "\t"; print "peak_set", "reads_in_peaks", "reads", "frip" }}
            NR == FNR {{ hits[$2] = $1; next }}
            {{ print $0, hits[$0] + 0, reads, (hits[$0] + 0) / reads }}' \
            "$tmp/hits" <(printf '%s\n' "${{labels[@]}}") \
            > {output:q}
        """


# ENCODE's JSD: deeptools fingerprints of the replicates and the first control, on
# BAMs without reads that touch an excluded region.
rule jsd:
    input:
        bams=lambda w: [f"{OUTDIR}/{w.exp}/filter/{w.exp}.{r}.bam" for r in EXPERIMENTS[w.exp]["libs"]],
        control=lambda w: (
            [f"{OUTDIR}/{EXPERIMENTS[w.exp]['control']}/filter/{EXPERIMENTS[w.exp]['control']}.ctl1.bam"]
            if EXPERIMENTS[w.exp]["control"]
            else []
        ),
        blacklist=BLACKLIST,
    output:
        qc=f"{OUTDIR}/{{exp}}/jsd/{{exp}}.tsv",
        plot=f"{OUTDIR}/{{exp}}/jsd/{{exp}}.svg",
        plot_png=f"{OUTDIR}/{{exp}}/jsd/{{exp}}.png",
    params:
        labels=lambda w: EXPERIMENTS[w.exp]["libs"] + ["ctl1"] * bool(EXPERIMENTS[w.exp]["control"]),
        mapq=lambda w: config["params"]["mapq_thresh"],
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/jsd/{{exp}}.txt",
    benchmark:
        f"{BENCHDIR}/jsd/{{exp}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        filtered=()
        for bam in {input.bams:q} {input.control:q}; do
            out="$tmp/${{#filtered[@]}}.bam"
            bedtools intersect -nonamecheck -v -abam "$bam" -b {input.blacklist:q} > "$out"
            samtools index "$out"
            filtered+=("$out")
        done
        controls=({input.control:q})
        jsd=()
        if (( ${{#controls[@]}} )); then
            jsd=(--JSDsample "${{filtered[-1]}}")
        fi
        plotFingerprint \
            -b "${{filtered[@]}}" \
            "${{jsd[@]}}" \
            --labels {params.labels:q} \
            --outQualityMetrics {output.qc:q} \
            --outRawCounts "$tmp/counts.tab" \
            --minMappingQuality {params.mapq} \
            --numberOfProcessors {threads} \
            --plotFile "$tmp/fingerprint.png"
        python {workflow.basedir:q}/scripts/plot_fingerprint.py \
            --counts "$tmp/counts.tab" \
            --title {wildcards.exp:q} \
            --plot {output.plot:q}
        """


# ENCODE's GC bias: Picard's metrics, replotted. ENCODE first strips read groups,
# which these BAMs never have.
rule gc_bias:
    input:
        bam=f"{OUTDIR}/{{exp}}/filter/{{exp}}.{{rep}}.bam",
        fasta=REFERENCE_FASTA,
    output:
        metrics=f"{OUTDIR}/{{exp}}/gc_bias/{{exp}}.{{rep}}.gc.txt",
        summary=f"{OUTDIR}/{{exp}}/gc_bias/{{exp}}.{{rep}}.gc_summary.txt",
        plot=f"{OUTDIR}/{{exp}}/gc_bias/{{exp}}.{{rep}}.svg",
        plot_png=f"{OUTDIR}/{{exp}}/gc_bias/{{exp}}.{{rep}}.png",
    params:
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/gc_bias/{{exp}}.{{rep}}.txt",
    benchmark:
        f"{BENCHDIR}/gc_bias/{{exp}}.{{rep}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        picard -Xmx4g CollectGcBiasMetrics \
            -R {input.fasta:q} \
            -I {input.bam:q} \
            -O {output.metrics:q} \
            -CHART "$tmp/gc.pdf" \
            -S {output.summary:q} \
            -ASSUME_SORTED false \
            -TMP_DIR "$tmp" \
            -VERBOSITY ERROR \
            -QUIET true
        python {workflow.basedir:q}/scripts/plot_gc.py --metrics {output.metrics:q} --plot {output.plot:q}
        """


def peak_qc_sets(exp):
    """Return the filtered peak sets whose widths and counts `peak_qc` reports.

    Parameters
    ----------
    exp : str
        ChIP experiment.

    Returns
    -------
    dict of str to str
        Path by label: each replicate and pool, then `<method>_opt` and
        `<method>_consv` unless `true_rep_only`.
    """
    sets = {b: peaks(exp, b, bfilt=True) for b in bases(exp)}
    if not config["params"]["true_rep_only"]:
        for method in methods(exp):
            for kind, short in (("optimal", "opt"), ("conservative", "consv")):
                sets[f"{method}_{short}"] = (
                    f"{OUTDIR}/{exp}/reproducibility/{exp}.{method}.{kind}.{peak_type(exp)}.gz"
                )
    return sets


# ENCODE's peak region size and peak count QC.
rule peak_qc:
    input:
        lambda w: list(peak_qc_sets(w.exp).values()),
    output:
        table=f"{OUTDIR}/{{exp}}/peak_qc/{{exp}}.tsv",
        plot=f"{OUTDIR}/{{exp}}/peak_qc/{{exp}}.svg",
        plot_png=f"{OUTDIR}/{{exp}}/peak_qc/{{exp}}.png",
    params:
        labels=lambda w: list(peak_qc_sets(w.exp)),
    log:
        f"{LOGDIR}/peak_qc/{{exp}}.txt",
    benchmark:
        f"{BENCHDIR}/peak_qc/{{exp}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/peak_qc.py \
            --peaks {input:q} \
            --labels {params.labels:q} \
            --title {wildcards.exp:q} \
            --table {output.table:q} \
            --plot {output.plot:q}
        """
