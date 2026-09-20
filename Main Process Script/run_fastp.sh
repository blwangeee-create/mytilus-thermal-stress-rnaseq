#!/bin/bash

RAW=/Volumes/Expansion/NCBI/PRJNA934294/raw_data
OUT=/Volumes/Expansion/NCBI/PRJNA934294/qc/fastp

for sample in SRR23438721 SRR23438722 SRR23438723 SRR23438724 SRR23438725 SRR23438726 SRR23438727 SRR23438728

do

echo "========== Processing ${sample} =========="

fastp \
-i ${RAW}/${sample}_1.fastq \
-I ${RAW}/${sample}_2.fastq \
-o ${OUT}/${sample}_clean_1.fastq \
-O ${OUT}/${sample}_clean_2.fastq \
-h ${OUT}/${sample}_fastp.html \
-j ${OUT}/${sample}_fastp.json \
-w 4

done
#!/bin/bash
set -euo pipefail

##########################################################
## Build Supplementary QC Tables
##########################################################

PROJECT=/Volumes/Expansion/NCBI/PRJNA934294

FASTP=${PROJECT}/qc/fastp
ALIGN=${PROJECT}/results/alignment
COUNT=${PROJECT}/results/counts/featurecounts_s0/counts_s0.txt.summary

OUT=${PROJECT}/results/qc_summary

mkdir -p ${OUT}

TABLE1=${OUT}/Supplementary_Table_S1_Final.tsv
TABLE2=${OUT}/Supplementary_Table_S2_Assignment.tsv

##########################################################
## Header
##########################################################

echo -e "SampleID\tRawReads\tCleanReads\tRetention(%)\tMeanReadLength(bp)\tQ20(%)\tQ30(%)\tGC(%)\tOverallMapping(%)\tAssignedReads\tAssigned(%)\tMultiMapping(%)\tUnmapped(%)" > ${TABLE1}

##########################################################
## Read featureCounts summary
##########################################################

declare -A Assigned
declare -A Multi
declare -A Unmapped
declare -A NoFeature
declare -A Ambiguous

header=($(head -1 ${COUNT}))

for ((i=1;i<${#header[@]};i++))
do
    sample=$(basename ${header[$i]} .sorted.bam)
    samples[$i]=$sample
done

while IFS=$'\t' read status values
do

    arr=($values)

    case ${status} in

    Assigned)
        for ((i=0;i<${#arr[@]};i++))
        do
            Assigned[${samples[$((i+1))]}]=${arr[$i]}
        done
    ;;

    Unassigned_MultiMapping)
        for ((i=0;i<${#arr[@]};i++))
        do
            Multi[${samples[$((i+1))]}]=${arr[$i]}
        done
    ;;

    Unassigned_Unmapped)
        for ((i=0;i<${#arr[@]};i++))
        do
            Unmapped[${samples[$((i+1))]}]=${arr[$i]}
        done
    ;;

    Unassigned_NoFeatures)
        for ((i=0;i<${#arr[@]};i++))
        do
            NoFeature[${samples[$((i+1))]}]=${arr[$i]}
        done
    ;;

    Unassigned_Ambiguity)
        for ((i=0;i<${#arr[@]};i++))
        do
            Ambiguous[${samples[$((i+1))]}]=${arr[$i]}
        done
    ;;

    esac

done < <(tail -n +2 ${COUNT})

##########################################################
## Table S1
##########################################################

for json in ${FASTP}/*_fastp.json
do

sample=$(basename ${json} _fastp.json)

raw=$(jq '.summary.before_filtering.total_reads' ${json})

clean=$(jq '.summary.after_filtering.total_reads' ${json})

q20=$(jq '.summary.after_filtering.q20_rate*100' ${json})

q30=$(jq '.summary.after_filtering.q30_rate*100' ${json})

gc=$(jq '.summary.after_filtering.gc_content*100' ${json})

len1=$(jq '.summary.after_filtering.read1_mean_length' ${json})

len2=$(jq '.summary.after_filtering.read2_mean_length' ${json})

length=$(( (len1+len2)/2 ))

retain=$(awk "BEGIN{printf \"%.2f\",100*${clean}/${raw}}")

map=$(grep "overall alignment rate" ${ALIGN}/${sample}_summary.txt | awk '{print $1}' | sed 's/%//')

assign=${Assigned[$sample]}

assignrate=$(awk "BEGIN{printf \"%.2f\",100*${assign}/${clean}}")

multi=${Multi[$sample]}

multirate=$(awk "BEGIN{printf \"%.2f\",100*${multi}/${clean}}")

unmap=${Unmapped[$sample]}

unmaprate=$(awk "BEGIN{printf \"%.2f\",100*${unmap}/${clean}}")

printf "%s\t%s\t%s\t%s\t%d\t%.2f\t%.2f\t%.2f\t%s\t%s\t%s\t%s\t%s\n" \
$sample \
$raw \
$clean \
$retain \
$length \
$q20 \
$q30 \
$gc \
$map \
$assign \
$assignrate \
$multirate \
$unmaprate >> ${TABLE1}

done

##########################################################
## Table S2
##########################################################

echo -e "SampleID\tAssigned\tUnmapped\tMultiMapping\tNoFeature\tAmbiguous" > ${TABLE2}

for s in $(printf "%s\n" ${!Assigned[@]} | sort)
do

echo -e "${s}\t${Assigned[$s]}\t${Unmapped[$s]}\t${Multi[$s]}\t${NoFeature[$s]}\t${Ambiguous[$s]}" >> ${TABLE2}

done

echo
echo "========================================"
echo "Finished!"
echo
echo ${TABLE1}
echo ${TABLE2}
echo "========================================"