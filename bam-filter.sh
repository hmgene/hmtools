
d=/mnt/vstor/SOM_GENE_BEG33/ChIP_seq/hg38/
s=TCCC-ST78_DMSO_H3K27ac_HiChIP_022924_CWRU_2
i=$d/$s/$s.bam

l=chr11:17,719,154-17,722,475


mamba activate corea
samtools view $i `echo $l | tr -d","`
