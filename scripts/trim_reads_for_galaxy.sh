#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# 0) Where your raw FASTQs live:
INPUT_DIR="/Users/tenanyehaglund/downloads/D142-445660387/BCL_Convert_03_26_2025_14_30_56-813607797/FASTQ/Original_FASTQ"

# 1) Where trimmed files will go (created alongside your repo):
OUTPUT_DIR="$PWD/trimmed_for_galaxy"
mkdir -p "$OUTPUT_DIR"

# 2) Trim parameters (same as you ran):
ADAPTER="AGATCGGAAGAGC"
QUALITY=25
MINLEN=50

# 3) Loop over every paired sample:
for R1 in "$INPUT_DIR"/CS_Lib*_R1_001.fastq.gz; do
  [[ -e "$R1" ]] || continue
  PREFIX=$(basename "$R1" _R1_001.fastq.gz)
  R2="$INPUT_DIR/${PREFIX}_R2_001.fastq.gz"
  if [[ ! -f "$R2" ]]; then
    echo "⚠️  No R2 for ${PREFIX}, skipping."
    continue
  fi

  OUT1="$OUTPUT_DIR/retrimmed_${PREFIX}_R1.fastq.gz"
  OUT2="$OUTPUT_DIR/retrimmed_${PREFIX}_R2.fastq.gz"

  echo "→ Trimming ${PREFIX}"
  cutadapt \
    -q ${QUALITY},${QUALITY} \
    -a ${ADAPTER} -A ${ADAPTER} \
    -m ${MINLEN} \
    -o "$OUT1" -p "$OUT2" \
    "$R1" "$R2"
done

echo "✅ All done. Trimmed files are in $OUTPUT_DIR"
