#!/bin/bash

# Set adapter sequence, quality threshold, and minimum length
ADAPTER="AGATCGGAAGAGC"
QUALITY=25
MIN_LENGTH=50

# Loop through all R1 files and process each paired-end sample
for R1 in CS_Lib*_R1_001.fastq.gz; do
    # Derive the matching R2 filename by replacing _R1_ with _R2_
    R2="${R1/_R1_/_R2_}"

    # Ensure R2 exists before running Cutadapt
    if [[ ! -f "$R2" ]]; then
        echo "Warning: Matching R2 file not found for $R1. Skipping..."
        continue
    fi

    # Define output filenames (prefix with "retrimmed_")
    OUT_R1="retrimmed_${R1}"
    OUT_R2="retrimmed_${R2}"

    echo "Processing sample: ${R1} and ${R2}"

    # Run Cutadapt with stricter quality trimming:
    cutadapt -q ${QUALITY},${QUALITY} \
        -a ${ADAPTER} -A ${ADAPTER} \
        -m ${MIN_LENGTH} \
        -o ${OUT_R1} -p ${OUT_R2} \
        ${R1} ${R2}

done
