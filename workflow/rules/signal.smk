# Signal tracks as bigWig, as the ENCODE ChIP-seq pipeline makes them.


# ENCODE's MACS2 signal tracks, one job per track. The bedGraphs are per million
# reads (--SPMR), so -S scales the p-value track back: treatment tagAlign lines over
# a million. Clamped to the contigs, sorted, and lines that start inside the one
# before dropped. slop clamps every interval past a contig end to its last base, so
# those lines tie: the sort's whole-line tiebreak and the awk pick the one kept.
rule signal:
    input:
        counts=lambda w: counts(w.exp, sample_ta(w.exp, w.base)),
        treat=f"{OUTDIR}/{{exp}}/macs3/{{exp}}.{{base}}.treat_pileup.bdg",
        control=f"{OUTDIR}/{{exp}}/macs3/{{exp}}.{{base}}.control_lambda.bdg",
        chromsizes=CHROMSIZES,
    output:
        f"{OUTDIR}/{{exp}}/signal/{{exp}}.{{base}}.{{track}}.bigWig",
    params:
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/signal/{{exp}}.{{base}}.{{track}}.txt",
    benchmark:
        f"{BENCHDIR}/signal/{{exp}}.{{base}}.{{track}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        if [[ {wildcards.track} == fc ]]; then
            macs3 bdgcmp -t {input.treat:q} -c {input.control:q} --outdir "$tmp" --o-prefix track -m FE
            bedgraph="$tmp/track_FE.bdg"
        else
            read -r lines _ < {input.counts:q}
            sval=$(awk -v n="$lines" 'BEGIN {{ printf "%.6f", n / 1000000 }}')
            macs3 bdgcmp -t {input.treat:q} -c {input.control:q} --outdir "$tmp" --o-prefix track -m ppois -S "$sval"
            bedgraph="$tmp/track_ppois.bdg"
        fi
        bedtools slop -i "$bedgraph" -g {input.chromsizes:q} -b 0 \
            | awk '$3 != -1' \
            | bedClip stdin {input.chromsizes:q} /dev/stdout \
            | sort -k1,1 -k2,2n -S 1G -T "$tmp" \
            | awk 'BEGIN {{ OFS = "\t" }} {{
                if (NR == 1 || prev_chr != $1 || prev_end <= $2) print
                prev_chr = $1; prev_end = $3
            }}' \
            > "$tmp/sorted.bdg"
        bedGraphToBigWig "$tmp/sorted.bdg" {input.chromsizes:q} {output:q}
        """


# ENCODE's count signal tracks: 5' read ends per strand. figwig counts every read
# in the filtered BAMs, the reads of the tagAligns before any subsampling. A pool
# reads every replicate's BAM.
rule count_signal:
    input:
        bams=lambda w: [
            f"{OUTDIR}/{w.exp}/filter/{w.exp}.{lib}.bam"
            for lib in (EXPERIMENTS[w.exp]["libs"] if w.base == "pooled" else [w.base])
        ],
        chromsizes=CHROMSIZES,
    output:
        positive=f"{OUTDIR}/{{exp}}/signal/{{exp}}.{{base}}.positive.bigWig",
        negative=f"{OUTDIR}/{{exp}}/signal/{{exp}}.{{base}}.negative.bigWig",
    params:
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/count_signal/{{exp}}.{{base}}.txt",
    benchmark:
        f"{BENCHDIR}/count_signal/{{exp}}.{{base}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        figwig bam2bw \
            -s {input.chromsizes:q} \
            -n "$tmp/counts" \
            -p {threads} \
            {input.bams:q}
        mv "$tmp/counts.+.bw" {output.positive:q}
        mv "$tmp/counts.-.bw" {output.negative:q}
        """
