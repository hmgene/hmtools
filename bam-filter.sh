#!/bin/bash
d=/mnt/vstor/SOM_GENE_BEG33/ChIP_seq/hg38/DATA
s=TCCC-ST78_DMSO_H3K27ac_HiChIP_022924_CWRU_2
i=$d/$s/$s.bam
o=$s.k3fil.bam

l=`echo chr11:17,719,154-17,722,475 | tr -d "," `
l="";

source activate corea
k=3
{
	samtools view -H $i
	samtools view $i $l | perl -ne '$_=~/NM:i:(\d+)/; if( $1 <= '$k'){ print $_;}' 
} | samtools view -b > $o
samtools index $o
