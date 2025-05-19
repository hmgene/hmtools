
cellranger-example(){
echo '
cellranger count --id=run_count_1kpbmcs \
   --fastqs=/mnt/home/user.name/yard/run_cellranger_count/pbmc_1k_v3_fastqs \
   --sample=pbmc_1k_v3 \
   --transcriptome=/mnt/home/user.name/yard/run_cellranger_count/refdata-gex-GRCh38-2020-A
'
}

cellranger-count(){
cellranger=$HMHOME/bigdata/cellranger-9.0.1/bin/cellranger
transcriptome=$HMHOME/bigdata/refdata-gex-GRCh38-2024-A
usage="$FUNCNAME <id> <fq> <sample> <transcriptome>"
if [ $# -lt 1 ];then echo "$usage";return;fi
	$cellranger count --id  $1 --fastqs=$2  --sample=$3 \
	   --transcriptome=$transcriptome --create-bam true
}

