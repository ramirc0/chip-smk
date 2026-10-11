#!/usr/bin/env bash
# Builds resources/fixture: the chr19 and chrM reads of two ENCODE ChIP-seq experiments
# and their controls, the chr19 + chrM genome, and the exclusion list cut to those contigs.
# Experiment reads come from the raw FASTQs; control reads, whose FASTQs are not on the
# mirror, from ENCODE's unfiltered BAMs, which hold every read.
#   YCR: ENCSR014YCR, ATF7 in GM12878, PE 100. Control ENCSR398JTO.
#   BGT: ENCSR000BGT, BATF in GM12878, SE 36. Control ENCSR000BGH, three replicates.
# Needs the Weng Lab ENCODE mirror. Run from the repository root:
#   pixi run --manifest-path workflow/envs/pixi.toml -e test bash .tests/integration/make_fixture.sh
set -euo pipefail

out=resources/fixture
encode=/data/projects/encode/data
fasta=/data/zusers/ramirezc/atac-smk/resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
blacklist=/data/zusers/ramirezc/atac-smk-refs/GRCh38/ENCFF356LFX.bed.gz
mkdir -p "$out"

samtools faidx "$fasta" chr19 chrM > "$out/genome.fa"
zcat "$blacklist" | awk '$1 == "chr19" || $1 == "chrM"' | pigz -n > "$out/ENCFF356LFX.bed.gz"

# Names of reads, or pairs by their first mate, that ENCODE placed on chr19, and every
# tenth on chrM.
names() {
    samtools view -@ 4 -F 0x980 "$encode/$1.bam" \
        | awk '$3 == "chr19" || ($3 == "chrM" && mito++ % 10 == 0) { print $1 }'
}

# name, BAM, then R1 and R2 FASTQs, or nothing to rebuild the reads from the BAM.
library() {
    local name=$1 bam=$2 r1=${3:-} r2=${4:-}
    local tmp
    tmp=$(mktemp -d -p "$out")
    names "$bam" > "$tmp/names"
    if [[ -n "$r1" ]]; then
        seqkit grep -j 4 -f "$tmp/names" "$encode/$r1.fastq.gz" -o "$out/${name}_R1.fastq.gz"
        [[ -z "$r2" ]] || seqkit grep -j 4 -f "$tmp/names" "$encode/$r2.fastq.gz" -o "$out/${name}_R2.fastq.gz"
    else
        samtools view -@ 4 -u -N "$tmp/names" "$encode/$bam.bam" \
            | samtools collate -u -O - "$tmp/collate" \
            | samtools fastq -n -1 "$tmp/R1.fastq" -2 "$tmp/R2.fastq" -0 "$tmp/R0.fastq" -s /dev/null -
        if [[ -s "$tmp/R1.fastq" ]]; then
            pigz -n < "$tmp/R1.fastq" > "$out/${name}_R1.fastq.gz"
            pigz -n < "$tmp/R2.fastq" > "$out/${name}_R2.fastq.gz"
        else
            pigz -n < "$tmp/R0.fastq" > "$out/${name}_R1.fastq.gz"
        fi
    fi
    rm -r "$tmp"
}
export -f names library
export out encode

parallel --colsep ' ' --halt now,fail=1 library {1} {2} {3} {4} <<'EOF2'
YCR_rep1 ENCSR014YCR/ENCFF767TBJ ENCSR014YCR/ENCFF044LTS ENCSR014YCR/ENCFF531UQG
YCR_rep2 ENCSR014YCR/ENCFF681XMB ENCSR014YCR/ENCFF317NXG ENCSR014YCR/ENCFF741BVE
YCR_ctl1 ENCSR398JTO/ENCFF611CPT
YCR_ctl2 ENCSR398JTO/ENCFF841YBE
BGT_rep1 ENCSR000BGT/ENCFF944ZNU ENCSR000BGT/ENCFF000NSW
BGT_rep2 ENCSR000BGT/ENCFF030DZL ENCSR000BGT/ENCFF000NSX
BGT_ctl1 ENCSR000BGH/ENCFF745UOD
BGT_ctl2 ENCSR000BGH/ENCFF379YGC
BGT_ctl3 ENCSR000BGH/ENCFF070THX
EOF2
