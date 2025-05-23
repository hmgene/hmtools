
bismark-merge-cov(){
usage="$FUNCNAME <input.txt: sample_name, file_name.cov.gz>

example usage:
for f in bigdata/methylseq/*/out/bismark/methylation_calls/methylation_coverage/*.deduplicated.bismark.cov.gz;do
    s=`echo $f | cut -d"/" -f 3`
    echo "$s $f"
done |  bismark-merge-cov - | gzip -c > merged_cov.gz

"
if [ $# -lt 1 ];then echo "$usage";return; fi
cat $1 | Rscript -e 'inp=read.table("stdin",header=F)
    library(data.table)
    d=NULL;
    for( i in 1:nrow(inp)){
        tt=fread(tt[i,2], col.names=c("chrom","start","end",paste(tt[i,1],c("perc","numC","numT"),sep="_")))
        if(is.null(d)){ d=tt;
        }else{ d=merge(d,tt,all=T); }
    }
    fwrite(d,"",sep="\t")
'
}


bismark-cov2bg-test(){
echo "
chr1	3000827	3000827	100	2	1
chr1	3000828	3000828	100	1	1
chr1	3001007	3001007	100	0	1
" | grep -v "^$" | hm bismark-cov2bg - 10
}

bismark2bg(){
echo "not implemtned"
#usage="$FUNCNAME <bam> <genome=mm10>";
#if [ $# -lt 2 ];then echo "$usage";return;fi
#
#bismark_methylation_extractor \
#  --bedGraph \
#  --cytosine_report \
#  --genome_folder $HMHOME/bigdata/mm10 \
#    bigdata/methylseq/bigdata/emseq-2025-04-11/2-year_1/out/bismark/deduplicated/2-year_1_val_1_bismark_bt2_pe.deduplicated.bam 
}
