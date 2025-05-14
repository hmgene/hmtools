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


pair2mcool ()
{
    c=$HMHOME/data/hg38.chrom.sizes;
    usage="$FUNCNAME <pair> [<binsz>=1k,5k,10k,20k,100k]";
    if [ $# -lt 1 ]; then
        echo "$usage";
        return;
    fi;
    local tmp=`mktemp -d`;
    local r=${2:-"1k,5k,10k,20k,100k"};
    r=`echo $r | tr [:lower:] [:upper:] | tr "," "\n" | numfmt --from=si | tr "\n" "," | sed "s/,$//"; `;
    cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 $c:${r%%,*} $1 $tmp/a.cool;
    cooler zoomify --balance -r$r $tmp/a.cool;
    cat $tmp/a.mcool
}
