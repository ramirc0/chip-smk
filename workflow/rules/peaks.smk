# Control choice, peak calls, naive overlap, IDR and the reproducible peak sets.

import math


# ENCODE's MACS2 call for ChIP-seq: tags extended to the fragment length, no shift.
MACS3_ARGS = [
    "-f", "BED",
    "-g", "hs",
    "-p", str(config["params"]["pval_thresh"]),
    "--nomodel",
    "--shift", "0",
    "--keep-dup", "all",
]

# ENCODE's MACS2 peak file: p-value order, Peak_N names, cap_num_peak rows, clipped
# to the contigs. The cap sits in awk because head would close the pipe early.
FORMAT_PEAKS = r"""
        LC_ALL=C sort -k8,8gr "$tmp/peaks_peaks.narrowPeak" \
            | awk -v cap={params.cap:q} 'BEGIN {{ OFS = "\t" }}
                NR <= cap {{
                    $4 = "Peak_" NR
                    if ($2 < 0) $2 = 0
                    if ($3 < 0) $3 = 0
                    if ($10 == -1) $10 = $2 + int(($3 - $2 + 1) / 2.0)
                    print
                }}' \
            > "$tmp/capped"
        bedClip -truncate "$tmp/capped" {input.chromsizes:q} "$tmp/clipped"
        pigz -n < "$tmp/clipped" > {output.peaks:q}
"""


# ENCODE's choose_ctl, from the tagAlign counts alone. A checkpoint: which control
# file each call reads, and whether it is a subsample, depends on read depths.
checkpoint choose_ctl:
    input:
        reps=lambda w: [counts(w.exp, effective(w.exp, r)) for r in EXPERIMENTS[w.exp]["libs"]],
        controls=lambda w: [
            counts(EXPERIMENTS[w.exp]["control"], effective(EXPERIMENTS[w.exp]["control"], c))
            for c in EXPERIMENTS[EXPERIMENTS[w.exp]["control"]]["libs"]
        ],
    output:
        f"{OUTDIR}/{{exp}}/choose_ctl/{{exp}}.tsv",
    params:
        flags=lambda w: choose_ctl_flags(w.exp),
    log:
        f"{LOGDIR}/choose_ctl/{{exp}}.txt",
    benchmark:
        f"{BENCHDIR}/choose_ctl/{{exp}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/choose_ctl.py \
            --reps {input.reps:q} \
            --controls {input.controls:q} \
            --output {output:q} \
            {params.flags:q}
        """


def choose_ctl_flags(exp):
    """Return the choose_ctl.py flag tokens from `config["params"]`.

    Parameters
    ----------
    exp : str
        ChIP experiment with a control.

    Returns
    -------
    list of str
        Flag tokens. Interpolate them with `:q`.
    """
    p = config["params"]
    flags = [
        "--ctl-depth-ratio", str(p["ctl_depth_ratio"]),
        "--ctl-depth-limit", str(p["ctl_depth_limit"]),
        "--exp-ctl-depth-ratio-limit", str(p["exp_ctl_depth_ratio_limit"]),
    ]
    if p["always_use_pooled_ctl"]:
        flags.append("--always-use-pooled-ctl")
    if layout(EXPERIMENTS[exp]["control"], "pooled") == "pe":
        flags.append("--paired-end-control")
    return flags


# ENCODE's SPP call, for tf experiments. run_spp.R would gunzip both tagAligns
# itself, single threaded.
rule spp:
    input:
        ta=lambda w: ta(w.exp, sample_ta(w.exp, w.sample)),
        control=lambda w: control_ta(w.exp, w.sample),
        fraglen=lambda w: fraglen(w.exp, w.sample),
        chromsizes=CHROMSIZES,
    output:
        peaks=f"{OUTDIR}/{{exp}}/spp/{{exp}}.{{sample}}.regionPeak.gz",
    params:
        cap=lambda w: config["params"]["cap_num_peak"]["spp"],
        fdr=lambda w: config["params"]["fdr_thresh"],
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/spp/{{exp}}.{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/spp/{{exp}}.{{sample}}.tsv"
    software:
        SPP_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        rapidgzip -d -c -P {threads} {input.ta:q} > "$tmp/chip.tagAlign"
        rapidgzip -d -c -P {threads} {input.control:q} > "$tmp/control.tagAlign"
        Rscript --max-ppsize=500000 "$(command -v run_spp.R)" \
            -c="$tmp/chip.tagAlign" \
            -i="$tmp/control.tagAlign" \
            -npeak={params.cap} \
            -odir="$tmp" \
            -speak="$(< {input.fraglen:q})" \
            -savr="$tmp/peaks.regionPeak" \
            -fdr={params.fdr} \
            -rf \
            -p={threads} \
            -tmpdir="$tmp"
        # Integer coordinates (SPP may write 1e+05), then clipped to the contigs.
        rapidgzip -d -c "$tmp/peaks.regionPeak.gz" \
            | awk 'BEGIN {{ OFS = "\t" }} {{ if ($2 < 0) $2 = 0; print $1, int($2), int($3), $4, $5, $6, $7, $8, $9, $10 }}' \
            > "$tmp/int"
        bedClip -truncate "$tmp/int" {input.chromsizes:q} "$tmp/clipped"
        pigz -n < "$tmp/clipped" > {output.peaks:q}
        """


# ENCODE runs this call twice per replicate and pool: for histone peaks, and with
# -B --SPMR for the signal tracks. One call gives both. tf experiments use only the
# bedGraphs.
rule macs3:
    input:
        ta=lambda w: ta(w.exp, sample_ta(w.exp, w.base)),
        control=lambda w: control_ta(w.exp, w.base),
        fraglen=lambda w: fraglen(w.exp, w.base),
        chromsizes=CHROMSIZES,
    output:
        peaks=f"{OUTDIR}/{{exp}}/macs3/{{exp}}.{{base}}.narrowPeak.gz",
        treat=temp(f"{OUTDIR}/{{exp}}/macs3/{{exp}}.{{base}}.treat_pileup.bdg"),
        control=temp(f"{OUTDIR}/{{exp}}/macs3/{{exp}}.{{base}}.control_lambda.bdg"),
    params:
        args=MACS3_ARGS,
        control=lambda w, input: ["-c", *input.control] if input.control else [],
        cap=lambda w: config["params"]["cap_num_peak"]["macs3"],
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/macs3/{{exp}}.{{base}}.txt",
    benchmark:
        f"{BENCHDIR}/macs3/{{exp}}.{{base}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        macs3 callpeak \
            -t {input.ta:q} \
            {params.control:q} \
            -n peaks \
            --outdir "$tmp" \
            {params.args:q} \
            --extsize "$(< {input.fraglen:q})" \
            -B \
            --SPMR
        mv "$tmp/peaks_treat_pileup.bdg" {output.treat:q}
        mv "$tmp/peaks_control_lambda.bdg" {output.control:q}
        """ + FORMAT_PEAKS


rule macs3_pseudorep:
    input:
        ta=lambda w: ta(w.exp, f"{w.base}-{w.pr}"),
        control=lambda w: control_ta(w.exp, w.base),
        fraglen=lambda w: fraglen(w.exp, w.base),
        chromsizes=CHROMSIZES,
    output:
        peaks=f"{OUTDIR}/{{exp}}/macs3/{{exp}}.{{base}}-{{pr}}.narrowPeak.gz",
    params:
        args=MACS3_ARGS,
        control=lambda w, input: ["-c", *input.control] if input.control else [],
        cap=lambda w: config["params"]["cap_num_peak"]["macs3"],
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/macs3_pseudorep/{{exp}}.{{base}}-{{pr}}.txt",
    benchmark:
        f"{BENCHDIR}/macs3_pseudorep/{{exp}}.{{base}}-{{pr}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        macs3 callpeak \
            -t {input.ta:q} \
            {params.control:q} \
            -n peaks \
            --outdir "$tmp" \
            {params.args:q} \
            --extsize "$(< {input.fraglen:q})"
        """ + FORMAT_PEAKS


# ENCODE's exclusion-list filter for every peak set: no contact with an excluded
# region, score capped at 1000, chr1-22, X and Y only (grep -P 'chr[\dXY]+\b').
rule blacklist_filter:
    input:
        peaks=f"{OUTDIR}/{{exp}}/{{peakdir}}/{{exp}}.{{unit}}.{{ptype}}.gz",
        blacklist=BLACKLIST,
    output:
        f"{OUTDIR}/{{exp}}/{{peakdir}}/{{exp}}.{{unit}}.bfilt.{{ptype}}.gz",
    log:
        f"{LOGDIR}/blacklist_filter/{{exp}}.{{peakdir}}.{{unit}}.{{ptype}}.txt",
    benchmark:
        f"{BENCHDIR}/blacklist_filter/{{exp}}.{{peakdir}}.{{unit}}.{{ptype}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bedtools intersect -nonamecheck -v -a {input.peaks:q} -b {input.blacklist:q} \
            | awk 'BEGIN {{ OFS = "\t" }} {{ if ($5 > 1000) $5 = 1000; print }}' \
            | awk '/chr[0-9XY]+([^A-Za-z0-9_]|$)/' \
            | pigz -n \
            > {output:q}
        """


def pair_inputs(exp, pair):
    """Return the raw peaks a pair compares: both members, then the oracle.

    Parameters
    ----------
    exp : str
        ChIP experiment.
    pair : str
        Pair from `pairs()`.

    Returns
    -------
    list of str
        Peak paths of the two members and of `oracle(pair)`.
    """
    return [peaks(exp, s) for s in (*pair.split("_vs_"), oracle(pair))]


# ENCODE's naive overlap: oracle peaks that overlap a peak in each member, by at
# least half of either peak.
rule overlap:
    input:
        lambda w: pair_inputs(w.exp, w.pair),
    output:
        temp(f"{OUTDIR}/{{exp}}/overlap/{{exp}}.{{pair}}.{{ptype}}.gz"),
    params:
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/overlap/{{exp}}.{{pair}}.{{ptype}}.txt",
    benchmark:
        f"{BENCHDIR}/overlap/{{exp}}.{{pair}}.{{ptype}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        set -- {input:q}
        # ENCODE decompresses first; bedtools mishandles two gzipped inputs.
        pigz -dc "$1" > "$tmp/peaks1"
        pigz -dc "$2" > "$tmp/peaks2"
        pigz -dc "$3" > "$tmp/oracle"
        half() {{
            awk 'BEGIN {{ FS = OFS = "\t" }} {{ s1 = $3 - $2; s2 = $13 - $12; if ($21 / s1 >= 0.5 || $21 / s2 >= 0.5) print }}' \
                | cut -f 1-10 \
                | sort \
                | uniq
        }}
        bedtools intersect -nonamecheck -wo -a "$tmp/oracle" -b "$tmp/peaks1" \
            | half \
            | bedtools intersect -nonamecheck -wo -a stdin -b "$tmp/peaks2" \
            | half \
            | pigz -n \
            > {output:q}
        """


# ENCODE's IDR, for tf experiments: SPP peaks ranked by signal value. idr needs 20
# merged peaks; with fewer, the outputs are empty instead of failing the run.
rule idr:
    input:
        peaks=lambda w: pair_inputs(w.exp, w.pair),
        chromsizes=CHROMSIZES,
    output:
        raw=f"{OUTDIR}/{{exp}}/idr/{{exp}}.{{pair}}.idr.txt.gz",
        peaks=temp(f"{OUTDIR}/{{exp}}/idr/{{exp}}.{{pair}}.regionPeak.gz"),
    params:
        threshold=lambda w: config["params"]["idr_thresh"],
        min_log10=lambda w: -math.log10(config["params"]["idr_thresh"]),
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/idr/{{exp}}.{{pair}}.txt",
    benchmark:
        f"{BENCHDIR}/idr/{{exp}}.{{pair}}.tsv"
    software:
        IDR_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        set -- {input.peaks:q}
        status=0
        idr \
            --samples "$1" "$2" \
            --peak-list "$3" \
            --input-file-type narrowPeak \
            --output-file "$tmp/idr" \
            --rank signal.value \
            --soft-idr-threshold {params.threshold:q} \
            --use-best-multisummit-IDR \
            2> "$tmp/stderr" \
            || status=$?
        cat "$tmp/stderr" >&2
        if (( status )); then
            grep -q 'at least 20 peaks post-merge' "$tmp/stderr"
            : > "$tmp/idr"
        fi
        bedClip -truncate "$tmp/idr" {input.chromsizes:q} "$tmp/clipped"
        pigz -n < "$tmp/clipped" > {output.raw:q}
        awk 'BEGIN {{ OFS = "\t" }} $12 >= {params.min_log10} {{
                if ($2 < 0) $2 = 0
                print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12
            }}' "$tmp/clipped" \
            | sort \
            | uniq \
            | sort -grk7,7 \
            | cut -f 1-10 \
            | pigz -n \
            > {output.peaks:q}
        """


rule plot_idr:
    input:
        f"{OUTDIR}/{{exp}}/idr/{{exp}}.{{pair}}.idr.txt.gz",
    output:
        svg=f"{OUTDIR}/{{exp}}/idr/{{exp}}.{{pair}}.svg",
        png=f"{OUTDIR}/{{exp}}/idr/{{exp}}.{{pair}}.png",
    params:
        threshold=lambda w: config["params"]["idr_thresh"],
    log:
        f"{LOGDIR}/plot_idr/{{exp}}.{{pair}}.txt",
    benchmark:
        f"{BENCHDIR}/plot_idr/{{exp}}.{{pair}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/plot_idr.py \
            --input {input:q} \
            --threshold {params.threshold:q} \
            --title "{wildcards.exp} {wildcards.pair}" \
            --output {output.svg:q}
        """


# ENCODE's reproducibility: the conservative set is the true pair with the most
# peaks, the optimal set that or the pooled pseudoreplicates' set, whichever is larger.
rule reproducibility:
    input:
        true_pairs=lambda w: [pair_peaks(w.exp, w.method, p) for p in pairs(w.exp) if "-" not in p],
        self_pseudo=lambda w: [
            pair_peaks(w.exp, w.method, f"{r}-pr1_vs_{r}-pr2") for r in EXPERIMENTS[w.exp]["libs"]
        ],
        pooled_pseudo=lambda w: [
            pair_peaks(w.exp, w.method, p) for p in pairs(w.exp) if p.startswith("pooled")
        ],
    output:
        optimal=f"{OUTDIR}/{{exp}}/reproducibility/{{exp}}.{{method}}.optimal.{{ptype}}.gz",
        conservative=f"{OUTDIR}/{{exp}}/reproducibility/{{exp}}.{{method}}.conservative.{{ptype}}.gz",
        table=f"{OUTDIR}/{{exp}}/reproducibility/{{exp}}.{{method}}.{{ptype}}.tsv",
    log:
        f"{LOGDIR}/reproducibility/{{exp}}.{{method}}.{{ptype}}.txt",
    benchmark:
        f"{BENCHDIR}/reproducibility/{{exp}}.{{method}}.{{ptype}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/peak_sets.py \
            --true-pairs {input.true_pairs:q} \
            --self-pseudo {input.self_pseudo:q} \
            --pooled-pseudo {input.pooled_pseudo:q} \
            --optimal {output.optimal:q} \
            --conservative {output.conservative:q} \
            --table {output.table:q}
        """
