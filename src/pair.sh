bam2pair () 
{ 
    c=$HMHOME/data/hg38.chrom.sizes;
    usage="$FUNCNAME <bam> [chrom=$c] > <pair>";
    if [ $# -lt 1 ]; then
        echo "$usage";
        return;
    fi;
    pairtools parse --chroms-path ${2:-$c} --drop-sam $1 | pairtools sort --nproc 8 - | pairtools dedup
}

