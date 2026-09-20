#!/bin/bash
set -euo pipefail

# ============================================================
# Build RNA-seq QC summary table
# ============================================================

PROJECT=/Volumes/Expansion/NCBI/PRJNA934294

FASTP_DIR=$PROJECT/qc/fastp
ALIGN_DIR=$PROJECT/results/alignment
COUNT_DIR=$PROJECT/results/counts/featurecounts_s0
OUT_DIR=$PROJECT/results/qc_summary

mkdir -p "$OUT_DIR"

OUTFILE=$OUT_DIR/Supplementary_Table_S1_QC.tsv

echo -e "SampleID\tRawReads\tCleanReads\tQ30(%)\tOverallMappingRate(%)\tAssignedReads\tAssignedRate(%)" > "$OUTFILE"

for json in "$FASTP_DIR"/*_fastp.json
do

sample=$(basename "$json" _fastp.json)

raw=$(jq '.summary.before_filtering.total_reads' "$json")

clean=$(jq '.summary.after_filtering.total_reads' "$json")

q30=$(jq -r '.summary.after_filtering.q30_rate*100' "$json")

map=$(grep "overall alignment rate" \
"$ALIGN_DIR/${sample}_summary.txt" | \
awk '{print $1}' | sed 's/%//')

assigned=$(awk -v s="$sample" '
NR==1{
for(i=2;i<=NF;i++)
if($i~s".sorted.bam") c=i
}
$1=="Assigned"{print $c}
' "$COUNT_DIR/counts_s0.txt.summary")

rate=$(awk -v a="$assigned" -v c="$clean" \
'BEGIN{printf "%.2f",a/c*100}')

printf "%s\t%s\t%s\t%.2f\t%s\t%s\t%s\n" \
"$sample" \
"$raw" \
"$clean" \
"$q30" \
"$map" \
"$assigned" \
"$rate" >> "$OUTFILE"

done

echo
echo "Finished!"
echo "$OUTFILE"
