#!/bin/bash

# yeast 3' UTR expression analysis
# (C) 2021 - Juan Caballero

# DEPS: fastp v020.1, hisat2 v2.2.1, stringtie v2.1.4, samtools v1.11, python v3.x

THREADS=8 # num of threads to use

# preprocess FastQ
fastp -w ${THREADS} -i SPADE.fastq.gz -o SPADE_trim.fastq.gz
mv fastp.html SPADE_trim_fastp_report.html
rm fastp.json

# build HISAT2 index for yeast
hisat2-build yeast_S88C.fasta yeast_S88C

# align reads to yeast genome
hisat2 --dta -x yeast_S88C -p ${THREADS} -U SPADE_trim.fastq.gz -S SPADE_align.sam

# generating sorted BAM
samtools view -bS SPADE_align.sam | samtools sort - > SPADE_align.bam
samtools index SPADE_align.bam
rm SPADE_align.sam

# detect 3'UTR
python detect3pUTR.py yeast_S88C.gff3 SPADE_align.bam yeast_3pUTR_predict.gff3

# do expression in 3' UTR
stringtie SPADE_align.bam -o SPADE_3pUTR_predict_annot.gtf -G yeast_3pUTR_predict.gff3 -p $THREADS -A 3pUTR_predict_expr.tsv

# rename chrom names, script expects [col num 0-based] [input] [output]
perl changeChrName.pl 2 3pUTR_predict_expr.tsv 3pUTR_predict_expr_chr.tsv