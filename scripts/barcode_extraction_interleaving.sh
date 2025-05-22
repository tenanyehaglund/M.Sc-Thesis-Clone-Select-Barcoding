#!/bin/bash
set -euo pipefail

# scCSlibrary3_Day8_2.sh
# Usage: bash ~/barcode_analysis/scCSlibrary3_Day8_2.sh

fwfastq="/Users/tenanyehaglund/downloads/D142-445660387/BCL_Convert_03_26_2025_14_30_56-813607797/FASTQ/Original_FASTQ/CS_Lib_Day8_2-TH-021325_S11_L001_R1_001.fastq.gz"
rvfastq="/Users/tenanyehaglund/downloads/D142-445660387/BCL_Convert_03_26_2025_14_30_56-813607797/FASTQ/Original_FASTQ/CS_Lib_Day8_2-TH-021325_S11_L001_R2_001.fastq.gz"
outdir="/Users/tenanyehaglund/barcode_analysis"

mkdir -p "$outdir"

# 1) Trim 20 bp before/after the 17 bp barcode AND interleave
cutadapt \
  -e 0.1 \
  --no-indels \
  -O 20 \
  -g "CTATTTCCGGTGAATTCCCG" \
  -G "AGCTCCTCGCCCTTGCTCAC" \
  -a "GTGAGCAAGGGCGAGGAGCT" \
  -A "CGGGAATTCACCGGAAATAG" \
  --discard-untrimmed \
  --interleaved \
  -o "$outdir/Day8_2_BClibmerge.fastq.gz" \
  "$fwfastq" \
  "$rvfastq"

echo "### Trimmed & interleaved → $outdir/Day8_2_BClibmerge.fastq.gz"

# 2) Decompress
gunzip -f "$outdir/Day8_2_BClibmerge.fastq.gz"
echo "### Decompressed → $outdir/Day8_2_BClibmerge.fastq"

# 3) Parse interleaved FASTQ to extract 17 nt from R1/R2 and count pairs
awk 'BEGIN {
    OFS = ","
}
NR % 8 == 2 {
    r1 = substr($0, 1, 17)
}
NR % 8 == 6 {
    r2 = substr($0, 1, 17)
    counts[r1 SUBSEP r2]++
}
END {
    print "barcode_R1","barcode_R2","count"
    for (key in counts) {
        split(key, a, SUBSEP)
        print a[1], a[2], counts[key]
    }
}' "$outdir/Day8_2_BClibmerge.fastq" \
  > "$outdir/Day8_2_barcodes_R1_R2_counts.csv"

echo "### Extracted → $outdir/Day8_2_barcodes_R1_R2_counts.csv"
