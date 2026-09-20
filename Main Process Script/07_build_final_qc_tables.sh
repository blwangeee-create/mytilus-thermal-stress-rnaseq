#!/bin/zsh
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
## 1. 用 awk 把 featureCounts 的 summary 转成易读表格
##########################################################

TMP_FC=${OUT}/tmp_fc_data.tsv

awk -F'\t' '
BEGIN { OFS="\t" }
NR==1 {
    # 记录样本数（列数减1）
    nsamples = NF - 1
    # 提取样本名
    for(i=2; i<=NF; i++) {
        n = split($i, path_parts, "/")
        base = path_parts[n]
        sub(/.sorted.bam$/, "", base)
        sample[i-1] = base
    }
}
NR>1 {
    status = $1
    # 把每一列的数值存入二维数组（用逗号做下标）
    for(i=2; i<=NF; i++) {
        value[status, i-1] = $i
    }
}
END {
    # 输出表头
    printf "SampleID\tAssigned\tUnmapped\tMultiMapping\tNoFeature\tAmbiguous\n"
    # 按样本顺序输出
    for(j=1; j<=nsamples; j++) {
        assigned   = value["Assigned", j]
        unmapped   = value["Unassigned_Unmapped", j]
        multi      = value["Unassigned_MultiMapping", j]
        nofeature  = value["Unassigned_NoFeatures", j]
        ambiguous  = value["Unassigned_Ambiguity", j]
        printf "%s\t%s\t%s\t%s\t%s\t%s\n", sample[j], assigned, unmapped, multi, nofeature, ambiguous
    }
}' ${COUNT} > ${TMP_FC}

# 把这个临时文件直接作为 Table S2
cp ${TMP_FC} ${TABLE2}

##########################################################
## 2. 生成 Table S1（整合 fastp、HISAT2、featureCounts）
##########################################################

echo -e "SampleID\tRawReads\tCleanReads\tRetention(%)\tMeanReadLength(bp)\tQ20(%)\tQ30(%)\tGC(%)\tOverallMapping(%)\tAssignedReads\tAssigned(%)\tMultiMapping(%)\tUnmapped(%)" > ${TABLE1}

for json in ${FASTP}/*_fastp.json
do
    sample=$(basename ${json} _fastp.json)

    # fastp 数据
    raw=$(jq '.summary.before_filtering.total_reads' ${json})
    clean=$(jq '.summary.after_filtering.total_reads' ${json})
    q20=$(jq '.summary.after_filtering.q20_rate*100' ${json})
    q30=$(jq '.summary.after_filtering.q30_rate*100' ${json})
    gc=$(jq '.summary.after_filtering.gc_content*100' ${json})
    len1=$(jq '.summary.after_filtering.read1_mean_length' ${json})
    len2=$(jq '.summary.after_filtering.read2_mean_length' ${json})
    length=$(( (len1+len2)/2 ))
    retain=$(awk "BEGIN{printf \"%.2f\",100*${clean}/${raw}}")

    # HISAT2 比对率
    map=$(grep "overall alignment rate" ${ALIGN}/${sample}_summary.txt | awk '{print $1}' | sed 's/%//')

    # 从临时文件中提取该样本的 featureCounts 统计
    fc_line=$(awk -F'\t' -v s="$sample" '$1 == s {print $0}' ${TMP_FC})
    assign=$(echo "$fc_line" | cut -f2)
    unmapped=$(echo "$fc_line" | cut -f3)
    multi=$(echo "$fc_line" | cut -f4)

    # 计算百分比
    assignrate=$(awk "BEGIN{printf \"%.2f\",100*${assign}/${clean}}")
    multirate=$(awk "BEGIN{printf \"%.2f\",100*${multi}/${clean}}")
    unmaprate=$(awk "BEGIN{printf \"%.2f\",100*${unmapped}/${clean}}")

    # 写入 Table S1
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

# 删除临时文件
rm ${TMP_FC}

echo
echo "========================================"
echo "Finished!"
echo
echo ${TABLE1}
echo ${TABLE2}
echo "========================================"
