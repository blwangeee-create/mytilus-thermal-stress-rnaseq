#!/bin/bash

# 严格定义核心路径（全部物理定位在外接硬盘）
INDEX="/Volumes/Expansion/NCBI/PRJNA934294/index/HISAT2_MCOR1.1/mcor1.1"
FASTQ_DIR="/Volumes/Expansion/NCBI/PRJNA934294/qc/fastp"
OUT_DIR="/Volumes/Expansion/NCBI/PRJNA934294/results/alignment"

mkdir -p ${OUT_DIR}

# 9个样本列表
SAMPLES=(SRR23438720 SRR23438721 SRR23438722 SRR23438723 SRR23438724 SRR23438725 SRR23438726 SRR23438727 SRR23438728)

echo "========================================="
echo "  开始 HISAT2 批量比对 + Samtools 实时压缩排序"
echo "========================================="

for SAMPLE in "${SAMPLES[@]}"; do
    echo "-----------------------------------------"
    echo "正在处理样本: ${SAMPLE} ..."
    
    # 工业级管道流：比对 -> 实时排序压缩
    hisat2 -p 4 \
           -x ${INDEX} \
           -1 ${FASTQ_DIR}/${SAMPLE}_clean_1.fastq \
           -2 ${FASTQ_DIR}/${SAMPLE}_clean_2.fastq \
           --summary-file ${OUT_DIR}/${SAMPLE}_summary.txt | \
    samtools sort -@ 2 -o ${OUT_DIR}/${SAMPLE}.sorted.bam -
    
    echo "样本 ${SAMPLE} 压缩排序完成！"
    echo "BAM文件：${OUT_DIR}/${SAMPLE}.sorted.bam"
    echo "比对报告：${OUT_DIR}/${SAMPLE}_summary.txt"
done

echo "========================================="
echo "  🎉 所有样本比对及流式压缩全部结束！"
echo "========================================="
