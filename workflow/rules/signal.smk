# Signal tracks as bigWig, as the ENCODE ChIP-seq pipeline makes them.


# ENCODE's MACS2 signal tracks, one job per track. The bedGraphs are per million
# reads (--SPMR), so -S scales the p-value track back: treatment tagAlign lines over
# a million. Clamped to the contigs, sorted, and lines that start inside the one
# before dropped.
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


# ENCODE's count signal tracks: 5' read ends per strand. One sort serves both.
rule count_signal:
    input:
        ta=lambda w: ta(w.exp, sample_ta(w.exp, w.base)),
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

        rapidgzip -d -c -P {threads} {input.ta:q} \
            | sort -k1,1 -k2,2n -S 2G --parallel {threads} -T "$tmp" \
            > "$tmp/sorted"
        for strand in positive negative; do
            sign=$([[ $strand == positive ]] && echo + || echo -)
            bedtools genomecov -5 -bg -strand "$sign" -g {input.chromsizes:q} -i "$tmp/sorted" > "$tmp/$strand.bdg"
        done
        bedGraphToBigWig "$tmp/positive.bdg" {input.chromsizes:q} {output.positive:q}
        bedGraphToBigWig "$tmp/negative.bdg" {input.chromsizes:q} {output.negative:q}
        """
