i=bwa-par 
. ../src/util.sh
echo "#!/bin/bash
hg38=/mnt/vstor/SOM_GENE_BEG33/ref_files/hg38/bwa_indices/GRChg38.d1.vd1.fa
ecoli=/mnt/vstor/SOM_GENE_BEG33/ref_files/Escherichia_coli_K_12_MG1655/Sequence/BWAIndex/genome.fa
dm3=/mnt/vstor/SOM_GENE_BEG33/ref_files/dm3/dm3.fa
mm10=/mnt/vstor/SOM_GENE_BEG33/ref_files/mm10/bwa_indices/mm10.fa
ercc=/mnt/vstor/SOM_GENE_BEG33/RNA_seq/hg19/ref/RSEM_ERCC/ERCC92.fa
genomes=( hg38 ecoli dm3 mm10 ercc)
export \${genomes[@]};
" > $i.sh
declare -f ck >> $i.sh
cat ../src/bwa.sh >> $i.sh
echo "bwa-parallel \$1 \$2 \$3" >> $i.sh

shc -U -v -r -f $i.sh -o $i


