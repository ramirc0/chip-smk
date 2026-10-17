# Reads to filtered BAMs and tagAligns, as the ENCODE ChIP-seq pipeline makes them.


# Shell functions shared by align and xcor_tagalign. r1 and r2 are bash arrays of
# FASTQs. bowtie2 takes lanes as comma lists and bwa reads them concatenated, so no
# lane is merged into a new file. --reorder keeps bowtie2's output in input order:
# records with equal coordinates keep it through the sort, and it seeds the shuffles.
ALIGN_FUNCTIONS = r"""
        lanes() {{ local IFS=,; echo "$*"; }}
        # One rapidgzip per lane: fed a BGZF lane without its EOF block and the next
        # lane as one stream, rapidgzip reports junk at the boundary and can hang.
        unzip_lanes() {{ local f; for f; do rapidgzip -d -c -P {threads} "$f"; done; }}
        # ENCODE's bwa check: drop every record of a read whose CIGAR, less deletions,
        # does not span its sequence. bwa writes a read's records together.
        drop_bad_cigar() {{
            awk 'BEGIN {{ FS = OFS = "\t" }}
                /^@/ {{ print; next }}
                $1 != name {{ if (!bad) printf "%s", group; name = $1; group = ""; bad = 0 }}
                {{ group = group $0 "\n" }}
                $6 != "*" {{
                    cigar = $6; gsub(/[0-9]+D/, "", cigar)
                    k = split(cigar, v, /[A-Z]/); s = 0
                    for (i = 1; i <= k; i++) s += v[i]
                    if (s != length($10)) bad = 1
                }}
                END {{ if (!bad) printf "%s", group }}'
        }}
        align() {{
            case {params.aligner}-$1 in
                bowtie2-se)
                    bowtie2 --mm --reorder --threads {threads} -x {params.index:q} \
                        -U "$(lanes "${{r1[@]}}")"
                    ;;
                bowtie2-pe)
                    bowtie2 -X 2000 --mm --reorder --threads {threads} -x {params.index:q} \
                        -1 "$(lanes "${{r1[@]}}")" -2 "$(lanes "${{r2[@]}}")"
                    ;;
                bwa-se)
                    bwa aln -q 5 -l 32 -k 2 -t {threads} {params.index:q} <(cat "${{r1[@]}}") > "$tmp/r1.sai"
                    bwa samse {params.index:q} "$tmp/r1.sai" <(cat "${{r1[@]}}")
                    ;;
                bwa-pe)
                    # ENCODE's read length: the longest of the first million R1 reads.
                    length=$(awk 'NR % 4 == 2 && length($0) > m {{ m = length($0) }} NR == 4000000 {{ exit }} END {{ print m + 0 }}' \
                        <(unzip_lanes "${{r1[@]}}"))
                    if (( {params.bwa_mem} && length >= {params.bwa_mem} )); then
                        bwa mem -M -t {threads} {params.index:q} <(cat "${{r1[@]}}") <(cat "${{r2[@]}}")
                    else
                        bwa aln -q 5 -l 32 -k 2 -t {threads} {params.index:q} <(cat "${{r1[@]}}") > "$tmp/r1.sai"
                        bwa aln -q 5 -l 32 -k 2 -t {threads} {params.index:q} <(cat "${{r2[@]}}") > "$tmp/r2.sai"
                        bwa sampe {params.index:q} "$tmp/r1.sai" "$tmp/r2.sai" <(cat "${{r1[@]}}") <(cat "${{r2[@]}}")
                    fi | drop_bad_cigar
                    ;;
            esac
        }}
"""


def align_params():
    """Return the params every aligning rule passes to `ALIGN_FUNCTIONS`.

    Returns
    -------
    dict of str to object
        `aligner`, `index` prefix and `bwa_mem`, the shortest R1 read length
        that switches paired-end bwa to `bwa mem` (0 = never).
    """
    return {
        "aligner": ALIGNER,
        "index": INDEX[ALIGNER],
        "bwa_mem": config["params"]["bwa_mem_min_read_len"],
    }


# Run-independent: built once under results/refs/, reused across runs.
rule chromsizes:
    input:
        REFERENCE_FASTA,
    output:
        CHROMSIZES,
    log:
        f"{LOGDIR}/chromsizes.txt",
    benchmark:
        f"{BENCHDIR}/chromsizes.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools faidx {input:q} -o - | cut -f 1,2 > {output:q}
        """


# conda's bowtie2-build lacks SAIS, so its size check alone picks a small .bt2
# index for GRCh38. --large-index keeps one index format for every genome.
rule bowtie2_build:
    input:
        REFERENCE_FASTA,
    output:
        INDEX_FILES["bowtie2"],
    params:
        prefix=INDEX["bowtie2"],
    log:
        f"{LOGDIR}/bowtie2_build.txt",
    benchmark:
        f"{BENCHDIR}/bowtie2_build.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bowtie2-build --large-index --threads {threads:q} {input:q} {params.prefix:q}
        """


rule bwa_index:
    input:
        REFERENCE_FASTA,
    output:
        INDEX_FILES["bwa"],
    params:
        prefix=INDEX["bwa"],
    log:
        f"{LOGDIR}/bwa_index.txt",
    benchmark:
        f"{BENCHDIR}/bwa_index.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bwa index -p {params.prefix:q} {input:q}
        """


# Report only, one lane at a time: every fastp trim and filter is off.
rule fastp_qc:
    input:
        r1=lambda w: lanes(w.exp, w.lib, "r1")[int(w.lane) - 1],
        r2=lambda w: lanes(w.exp, w.lib, "r2")[int(w.lane) - 1:int(w.lane)],
    output:
        html=f"{OUTDIR}/{{exp}}/fastp/{{exp}}.{{lib}}.{{lane}}.html",
        json=f"{OUTDIR}/{{exp}}/fastp/{{exp}}.{{lib}}.{{lane}}.json",
    params:
        r2=lambda w, input: ["-I", input.r2[0]] if input.r2 else [],
    log:
        f"{LOGDIR}/fastp_qc/{{exp}}.{{lib}}.{{lane}}.txt",
    benchmark:
        f"{BENCHDIR}/fastp_qc/{{exp}}.{{lib}}.{{lane}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        fastp \
            -i {input.r1:q} \
            {params.r2:q} \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --disable_length_filtering \
            --disable_trim_poly_g \
            -h {output.html:q} \
            -j {output.json:q} \
            -w {threads}
        """


# ENCODE crops with Trimmomatic MINLEN:<crop - tol> CROP:<crop>. cutadapt truncates
# (--length) before it filters (-m), so it keeps the same reads with the same bases.
# Pairs go when either mate does, as Trimmomatic's paired output. cutadapt cannot
# sniff the compression of a pipe, so it reads decompressed lanes.
rule align:
    input:
        r1=lambda w: lanes(w.exp, w.lib, "r1"),
        r2=lambda w: lanes(w.exp, w.lib, "r2"),
        index=INDEX_FILES[ALIGNER],
    output:
        bam=f"{OUTDIR}/{{exp}}/align/{{exp}}.{{lib}}.bam",
        flagstat=f"{OUTDIR}/{{exp}}/align/{{exp}}.{{lib}}.flagstat.json",
    params:
        **align_params(),
        layout=lambda w: LIBRARIES[w.exp, w.lib]["layout"],
        crop=lambda w: EXPERIMENTS[w.exp]["crop"],
        min_length=lambda w: max(EXPERIMENTS[w.exp]["crop"] - abs(config["params"]["crop_length_tol"]), 0),
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/align/{{exp}}.{{lib}}.txt",
    benchmark:
        f"{BENCHDIR}/align/{{exp}}.{{lib}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT
        """ + ALIGN_FUNCTIONS + r"""
        r1=({input.r1:q})
        r2=({input.r2:q})
        if (( {params.crop} )); then
            if [[ {params.layout} == pe ]]; then
                cutadapt -j {threads} -Z --length {params.crop} -m {params.min_length} \
                    -o "$tmp/R1.fastq.gz" -p "$tmp/R2.fastq.gz" \
                    <(unzip_lanes "${{r1[@]}}") <(unzip_lanes "${{r2[@]}}")
                r2=("$tmp/R2.fastq.gz")
            else
                cutadapt -j {threads} -Z --length {params.crop} -m {params.min_length} \
                    -o "$tmp/R1.fastq.gz" <(unzip_lanes "${{r1[@]}}")
            fi
            r1=("$tmp/R1.fastq.gz")
        fi

        align {params.layout} \
            | samtools sort -@ {threads} -m 1536M -T "$tmp/sort" -o {output.bam:q} -
        samtools flagstat -@ {threads} --output-fmt json {output.bam:q} > {output.flagstat:q}
        """


# ENCODE's unique-read filter: -F 1804 drops unmapped, secondary, QC-fail and
# duplicate records, -q the low MAPQ ones. Paired-end reads must stay properly
# paired (-f 2) after fixmate. filt holds what duplicate marking sees, name sorted
# for paired-end library complexity. chrM stays, as ENCODE's filter_chrs is empty.
rule filter:
    input:
        f"{OUTDIR}/{{exp}}/align/{{exp}}.{{lib}}.bam",
    output:
        bam=f"{OUTDIR}/{{exp}}/filter/{{exp}}.{{lib}}.bam",
        bai=f"{OUTDIR}/{{exp}}/filter/{{exp}}.{{lib}}.bam.bai",
        flagstat=f"{OUTDIR}/{{exp}}/filter/{{exp}}.{{lib}}.flagstat.json",
        dup=f"{OUTDIR}/{{exp}}/filter/{{exp}}.{{lib}}.dup.flagstat.json",
        filt=temp(f"{OUTDIR}/{{exp}}/filter/{{exp}}.{{lib}}.filt.bam"),
    params:
        layout=lambda w: LIBRARIES[w.exp, w.lib]["layout"],
        mapq=lambda w: config["params"]["mapq_thresh"],
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/filter/{{exp}}.{{lib}}.txt",
    benchmark:
        f"{BENCHDIR}/filter/{{exp}}.{{lib}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        if [[ {params.layout} == pe ]]; then
            # A spill to disk doubles the name sort's time.
            samtools view -@ {threads} -u -F 1804 -f 2 -q {params.mapq} {input:q} \
                | samtools sort -@ {threads} -m 2560M -n -u -T "$tmp/nmsrt" - \
                | samtools fixmate -@ {threads} -r -u - - \
                | samtools view -@ {threads} -1 -F 1804 -f 2 -o {output.filt:q} -
            samtools sort -@ {threads} -m 1536M -l 1 -T "$tmp/sort" -o "$tmp/filt.bam" {output.filt:q}
            filt="$tmp/filt.bam"
            proper=(-f 2)
        else
            # The aligned BAM is already in coordinate order.
            samtools view -@ {threads} -1 -F 1804 -q {params.mapq} -o {output.filt:q} {input:q}
            filt={output.filt:q}
            proper=()
        fi
        # ENCODE's sizes. They change which copy of a duplicate survives.
        sambamba markdup \
            -t {threads} \
            --tmpdir "$tmp" \
            --hash-table-size=17592186044416 \
            --overflow-list-size=20000000 \
            --io-buffer-size=256 \
            "$filt" \
            "$tmp/dupmark.bam"
        samtools flagstat -@ {threads} --output-fmt json "$tmp/dupmark.bam" > {output.dup:q}
        samtools view -@ {threads} -b -F 1804 "${{proper[@]}}" -o {output.bam:q} "$tmp/dupmark.bam"
        samtools index -@ {threads} {output.bam:q}
        samtools flagstat -@ {threads} --output-fmt json {output.bam:q} > {output.flagstat:q}
        """


# ENCODE's bam2ta: single-end reads in coordinate order; paired-end reads in name
# order, one line per mate, mate 1 first.
rule tagalign:
    input:
        f"{OUTDIR}/{{exp}}/filter/{{exp}}.{{lib}}.bam",
    output:
        ta=f"{OUTDIR}/{{exp}}/tagalign/{{exp}}.{{lib}}.tagAlign.gz",
        counts=f"{OUTDIR}/{{exp}}/tagalign/{{exp}}.{{lib}}.tagAlign.counts",
    params:
        layout=lambda w: LIBRARIES[w.exp, w.lib]["layout"],
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/tagalign/{{exp}}.{{lib}}.txt",
    benchmark:
        f"{BENCHDIR}/tagalign/{{exp}}.{{lib}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT
        """ + TA_FUNCTIONS + r"""
        if [[ {params.layout} == pe ]]; then
            samtools sort -@ {threads} -m 1536M -n -u -T "$tmp/nmsrt" {input:q} \
                | bedtools bamtobed -bedpe -mate1 -i stdin \
                | awk 'BEGIN {{ OFS = "\t" }} {{ print $1, $2, $3, "N", 1000, $9; print $4, $5, $6, "N", 1000, $10 }}'
        else
            bedtools bamtobed -i {input:q} \
                | bed_to_ta
        fi \
            | count {output.counts:q} \
            | pigz -n -p {threads} \
            > {output.ta:q}
        """


# ENCODE's subsample: a seeded shuffle of the whole file, keeping the first n reads
# (n / 2 pairs). The seed is the tagAlign's uncompressed size. shuf reads a pipe,
# as in ENCODE: on a small regular file, GNU shuf -n takes another code path.
rule subsample_ta:
    input:
        ta=f"{OUTDIR}/{{exp}}/tagalign/{{exp}}.{{name}}.tagAlign.gz",
        counts=f"{OUTDIR}/{{exp}}/tagalign/{{exp}}.{{name}}.tagAlign.counts",
    output:
        ta=f"{OUTDIR}/{{exp}}/tagalign/{{exp}}.{{name}}.sub{{n}}.tagAlign.gz",
        counts=f"{OUTDIR}/{{exp}}/tagalign/{{exp}}.{{name}}.sub{{n}}.tagAlign.counts",
    params:
        layout=lambda w: layout(w.exp, w.name),
    log:
        f"{LOGDIR}/subsample_ta/{{exp}}.{{name}}.sub{{n}}.txt",
    benchmark:
        f"{BENCHDIR}/subsample_ta/{{exp}}.{{name}}.sub{{n}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        """ + TA_FUNCTIONS + r"""
        read -r lines bytes < {input.counts:q}
        if [[ {params.layout} == pe ]]; then
            rapidgzip -d -c -P {threads} {input.ta:q} \
                | paste - - \
                | shuf -n $(( {wildcards.n} / 2 )) --random-source=<(keystream "$bytes") \
                | unpair
        else
            rapidgzip -d -c -P {threads} {input.ta:q} \
                | shuf -n {wildcards.n} --random-source=<(keystream "$bytes")
        fi \
            | count {output.counts:q} \
            | pigz -n -p {threads} \
            > {output.ta:q}
        """


# ENCODE's pool_ta: the members' reads, concatenated. gzip members concatenate.
rule pool_ta:
    input:
        tas=lambda w: [ta(w.exp, name) for name in pool_members(w.exp, w.pool)],
        counts=lambda w: [counts(w.exp, name) for name in pool_members(w.exp, w.pool)],
    output:
        ta=f"{OUTDIR}/{{exp}}/tagalign/{{exp}}.{{pool}}.tagAlign.gz",
        counts=f"{OUTDIR}/{{exp}}/tagalign/{{exp}}.{{pool}}.tagAlign.counts",
    log:
        f"{LOGDIR}/pool_ta/{{exp}}.{{pool}}.txt",
    benchmark:
        f"{BENCHDIR}/pool_ta/{{exp}}.{{pool}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        cat {input.tas:q} > {output.ta:q}
        awk '{{ n += $1; b += $2 }} END {{ printf "%d\t%d\n", n, b }}' {input.counts:q} > {output.counts:q}
        """


# ENCODE's spr: reads (whole pairs) shuffled with the seeded keystream, then halved.
rule pseudoreps:
    input:
        ta=lambda w: ta(w.exp, effective(w.exp, w.rep)),
        counts=lambda w: counts(w.exp, effective(w.exp, w.rep)),
    output:
        pr1=f"{OUTDIR}/{{exp}}/tagalign/{{exp}}.{{rep}}-pr1.tagAlign.gz",
        pr2=f"{OUTDIR}/{{exp}}/tagalign/{{exp}}.{{rep}}-pr2.tagAlign.gz",
        counts1=f"{OUTDIR}/{{exp}}/tagalign/{{exp}}.{{rep}}-pr1.tagAlign.counts",
        counts2=f"{OUTDIR}/{{exp}}/tagalign/{{exp}}.{{rep}}-pr2.tagAlign.counts",
    params:
        layout=lambda w: LIBRARIES[w.exp, w.rep]["layout"],
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/pseudoreps/{{exp}}.{{rep}}.txt",
    benchmark:
        f"{BENCHDIR}/pseudoreps/{{exp}}.{{rep}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT
        """ + TA_FUNCTIONS + r"""
        read -r lines bytes < {input.counts:q}
        if [[ {params.layout} == pe ]]; then
            rapidgzip -d -c -P {threads} {input.ta:q} \
                | paste - - \
                | shuf --random-source=<(keystream "$bytes") \
                | split -d -l $(( (lines / 2 + 1) / 2 )) - "$tmp/pr"
            half() {{ unpair < "$1"; }}
        else
            rapidgzip -d -c -P {threads} {input.ta:q} \
                | shuf --random-source=<(keystream "$bytes") \
                | split -d -l $(( (lines + 1) / 2 )) - "$tmp/pr"
            half() {{ cat "$1"; }}
        fi
        half "$tmp/pr00" | count {output.counts1:q} | pigz -n -p {threads} > {output.pr1:q}
        half "$tmp/pr01" | count {output.counts2:q} | pigz -n -p {threads} > {output.pr2:q}
        """
