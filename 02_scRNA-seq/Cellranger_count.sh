#!/bin/bash
FASTQ_DIR="/path/to/your/fastq/directory"
TRANSCRIPTOME_DIR="/path/to/your/reference/refdata-gex-GRCh38-2024-A"
OUTPUT_DIR="/path/to/your/output/directory"
SAMPLES=("SC-1" "SC-2" "SC-3" "SC-4" "SC-5")
cd ${OUTPUT_DIR} || exit
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Starting Cell Ranger count for sample: ${SAMPLE}"
    cellranger count \
        --id=${SAMPLE} \
        --fastqs=${FASTQ_DIR}/${SAMPLE} \
        --sample=${SAMPLE} \
        --transcriptome=${TRANSCRIPTOME_DIR} \
        --localcores=16 \
        --localmem=64 \
        --create-bam=true
    echo "Finished processing ${SAMPLE}."
done
