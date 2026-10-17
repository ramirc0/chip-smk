# Fragment length by strand cross-correlation, as ENCODE estimates it.


# ENCODE's xcor reads: R1 alone, trimmed to xcor_trim_bp, aligned single end, then
# MAPQ-filtered without duplicate removal. trimfastq.py trims every read or, when one
# read is shorter, none; it also writes "." bases as N. Never cropped. ENCODE also
# marks duplicates and computes QC on these reads, then discards both.
rule xcor_tagalign:
    input:
        r1=lambda w: lanes(w.exp, w.rep, "r1"),
        index=INDEX_FILES[ALIGNER],
    output:
        ta=temp(f"{OUTDIR}/{{exp}}/xcor/{{exp}}.{{rep}}.tagAlign.gz"),
        counts=temp(f"{OUTDIR}/{{exp}}/xcor/{{exp}}.{{rep}}.tagAlign.counts"),
    params:
        **align_params(),
        trim=lambda w: config["params"]["xcor_trim_bp"],
        mapq=lambda w: config["params"]["mapq_thresh"],
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/xcor_tagalign/{{exp}}.{{rep}}.txt",
    benchmark:
        f"{BENCHDIR}/xcor_tagalign/{{exp}}.{{rep}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT
        """ + ALIGN_FUNCTIONS + TA_FUNCTIONS + r"""
        r1=({input.r1:q})
        short=$(awk -v n={params.trim} 'NR % 4 == 2 && length($0) < n {{ print 1; exit }}' \
            <(unzip_lanes "${{r1[@]}}"))
        if [[ -z "$short" ]]; then
            unzip_lanes "${{r1[@]}}" \
                | awk -v n={params.trim} '
                    NR % 4 == 2 {{ $0 = substr($0, 1, n); gsub(/\./, "N") }}
                    NR % 4 == 0 {{ $0 = substr($0, 1, n) }}
                    {{ print }}' \
                | pigz -n -p {threads} \
                > "$tmp/R1.fastq.gz"
            r1=("$tmp/R1.fastq.gz")
        else
            echo "A read is shorter than {params.trim} bp: no trimming"
        fi

        align se \
            | samtools view -@ {threads} -u -F 1804 -q {params.mapq} - \
            | samtools sort -@ {threads} -m 1536M -u -T "$tmp/sort" - \
            | bedtools bamtobed -i stdin \
            | bed_to_ta \
            | count {output.counts:q} \
            | pigz -n -p {threads} \
            > {output.ta:q}
        """


# ENCODE's xcor: drop chrM, subsample with the keystream seeded by the whole
# tagAlign's size, then run_spp.R. The exclusion range ends at the read length,
# averaged over the first 100 subsampled reads, plus 10. -savd keeps the profile.
rule xcor:
    input:
        ta=f"{OUTDIR}/{{exp}}/xcor/{{exp}}.{{rep}}.tagAlign.gz",
        counts=f"{OUTDIR}/{{exp}}/xcor/{{exp}}.{{rep}}.tagAlign.counts",
    output:
        qc=f"{OUTDIR}/{{exp}}/xcor/{{exp}}.{{rep}}.cc.qc",
        fraglen=f"{OUTDIR}/{{exp}}/xcor/{{exp}}.{{rep}}.fraglen",
        profile=f"{OUTDIR}/{{exp}}/xcor/{{exp}}.{{rep}}.cc.tsv",
    params:
        mito=lambda w: config["mito_chr_name"],
        subsample=lambda w: config["params"]["xcor_subsample_reads"],
        min_max=lambda w: 50 if EXPERIMENTS[w.exp]["type"] == "tf" else 100,
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/xcor/{{exp}}.{{rep}}.txt",
    benchmark:
        f"{BENCHDIR}/xcor/{{exp}}.{{rep}}.tsv"
    software:
        SPP_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT
        """ + TA_FUNCTIONS + r"""
        read -r lines bytes < {input.counts:q}
        rapidgzip -d -c -P {threads} {input.ta:q} \
            | awk -v m={params.mito:q} '$1 != m' \
            | shuf -n {params.subsample} --random-source=<(keystream "$bytes") \
            > "$tmp/sub.tagAlign"
        length=$(head -n 100 "$tmp/sub.tagAlign" \
            | awk '{{ d = $3 - $2; s += d < 0 ? -d : d }} END {{ print int(s / NR) }}')
        end=$(( length + 10 > {params.min_max} ? length + 10 : {params.min_max} ))
        Rscript --max-ppsize=500000 "$(command -v run_spp.R)" \
            -rf \
            -c="$tmp/sub.tagAlign" \
            -p={threads} \
            -filtchr={params.mito:q} \
            -savd="$tmp/cc.RData" \
            -out="$tmp/cc.qc" \
            -x=-500:"$end" \
            -tmpdir="$tmp"
        # The first estimate of each comma-separated list, as ENCODE keeps.
        sed -r 's/,[^\t]+//g' "$tmp/cc.qc" > {output.qc:q}
        cut -f 3 {output.qc:q} > {output.fraglen:q}
        Rscript -e 'a <- commandArgs(TRUE); load(a[1]); write.table(crosscorr$cross.correlation, a[2], sep = "\t", quote = FALSE, row.names = FALSE)' \
            "$tmp/cc.RData" {output.profile:q}
        """


rule plot_xcor:
    input:
        qc=f"{OUTDIR}/{{exp}}/xcor/{{exp}}.{{rep}}.cc.qc",
        profile=f"{OUTDIR}/{{exp}}/xcor/{{exp}}.{{rep}}.cc.tsv",
    output:
        svg=f"{OUTDIR}/{{exp}}/xcor/{{exp}}.{{rep}}.svg",
        png=f"{OUTDIR}/{{exp}}/xcor/{{exp}}.{{rep}}.png",
    params:
        title=lambda w: f"{w.exp} {w.rep}",
    log:
        f"{LOGDIR}/plot_xcor/{{exp}}.{{rep}}.txt",
    benchmark:
        f"{BENCHDIR}/plot_xcor/{{exp}}.{{rep}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/plot_xcor.py \
            --qc {input.qc:q} \
            --profile {input.profile:q} \
            --title {params.title:q} \
            --output {output.svg:q}
        """


# The pooled replicate's fragment length: the replicates' mean, rounded half up as
# ENCODE's container rounds it (Python 2).
rule fraglen_mean:
    input:
        lambda w: [fraglen(w.exp, rep) for rep in EXPERIMENTS[w.exp]["libs"]],
    output:
        f"{OUTDIR}/{{exp}}/xcor/{{exp}}.pooled.fraglen",
    log:
        f"{LOGDIR}/fraglen_mean/{{exp}}.txt",
    benchmark:
        f"{BENCHDIR}/fraglen_mean/{{exp}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        cat {input:q} | awk '{{ s += $1 }} END {{ print int(s / NR + 0.5) }}' > {output:q}
        """
