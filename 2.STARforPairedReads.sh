#! /bin/bash

########## Get chimeric sites from paired-end reads

#### mapping for paired-end reads
starres=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/7.star
fileName=WM_SGWTS_293T_S_S1_L001
STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName} --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R1_001.filter3.fastq ${trimdir}/${fileName}_R2_001.filter3.fastq --runThreadN 10 --chimOutType WithinBAM --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1 
##
#samtools view -h ${starres}/${fileName}Aligned.sortedByCoord.out.bam | grep -E '(^@.*|ch:A:1)' | samtools view -Sb - > ${starres}/${fileName}.chimOut.bam
samtools view -h ${starres}/${fileName}Aligned.sortedByCoord.out.bam | awk 'BEGIN {FS="\t"} {if ($0 ~ /^@/ || ($0 ~ /ch:A:1/ && $0 ~ /HI:i:1/)) print $0}' | samtools view -bS - > ${starres}/${fileName}.chimOut.bam
bedtools bamtobed -i ${starres}/${fileName}.chimOut.bam -tag NM -cigar > ${starres}/${fileName}.chimOut.bed

#### extract bait and translocation sites
cd /mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/run
R CMD BATCH ./starforchimeric.sample.r


