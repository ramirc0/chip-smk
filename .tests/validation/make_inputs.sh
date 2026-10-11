#!/usr/bin/env bash
# Rebuilds the FASTQs of the validation libraries whose raw reads are not on the
# Weng Lab ENCODE mirror, from ENCODE's unfiltered BAMs. Those BAMs hold every
# read, unmapped ones included, after ENCODE's crop, so the sheet leaves their
# crop_length empty. Their reads come out in coordinate order, not ENCODE's, so the
# seeded xcor subsample and pseudoreplicates differ a little from ENCODE's draws.
# Run from the repository root:
#   pixi run --manifest-path workflow/envs/pixi.toml -e test bash .tests/validation/make_inputs.sh
set -euo pipefail

out=/data/zusers/ramirezc/chip-smk-refs/validation
encode=/data/projects/encode/data
mkdir -p "$out"

# layout (se or pe), then <experiment>/<unfiltered BAM accession>.
rebuild() {
    local layout=$1 bam=$2 tmp
    local name=${bam#*/}
    tmp=$(mktemp -d -p "$out")
    if [[ $layout == pe ]]; then
        # samtools fastq pairs mates only when they are adjacent.
        samtools collate -@ 4 -u -O "$encode/$bam.bam" "$tmp/collate" \
            | samtools fastq -@ 4 -n -1 "$tmp/R1.fastq.gz" -2 "$tmp/R2.fastq.gz" -0 /dev/null -s /dev/null -
        mv "$tmp/R2.fastq.gz" "$out/${name}_R2.fastq.gz"
    else
        samtools fastq -@ 4 -n -0 "$tmp/R1.fastq.gz" "$encode/$bam.bam"
    fi
    mv "$tmp/R1.fastq.gz" "$out/${name}_R1.fastq.gz"
    rm -r "$tmp"
}
export -f rebuild
export out encode

parallel --colsep ' ' --halt now,fail=1 -j 6 rebuild {1} {2} <<'LIST'
pe ENCSR398JTO/ENCFF611CPT
pe ENCSR398JTO/ENCFF841YBE
se ENCSR000BGH/ENCFF745UOD
se ENCSR000BGH/ENCFF379YGC
se ENCSR000BGH/ENCFF070THX
pe ENCSR736ZEG/ENCFF587OER
pe ENCSR736ZEG/ENCFF116ZVL
pe ENCSR982AYN/ENCFF216JKP
pe ENCSR982AYN/ENCFF845ZAN
se ENCSR243BPI/ENCFF941QLJ
se ENCSR243BPI/ENCFF902AQN
se ENCSR346XJA/ENCFF093EYW
LIST
