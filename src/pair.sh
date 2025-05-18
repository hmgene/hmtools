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


pair-intersect ()
{
    usage="$FUNCNAME <pair> <bed> <binsz> [options]"
    if [ $# -lt 3 ];then echo "$usage";return;fi;
    #columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
    #SRR30820835.49207559	chr1	10158	chr1	10187	+	-	UU
    local tmp=`mktemp -d `;
    hm mycat $1 > $tmp/a
    grep -v "^#" $tmp/a | perl -ne 'my $B='$3';chomp;my @d=split/\t/,$_;
        my ($b1,$b2)=map { int($_/$B)*$B } ($d[2],$d[4]);
        print join("\t",$d[1],$b1,$b1+$B,$d[3],$b2,$b2+$B,join("@",@d)),"\n";
    ' | bedtools pairtobed -a stdin -b $2 ${@:4} | cut -f 7 | tr "@" "\t" > $tmp/b
    grep "^#" $tmp/a
    cat $tmp/b
}



pair2mcool ()
{
    c=$HMHOME/data/hg38.chrom.sizes;
    usage="$FUNCNAME <pair> [<binsz>=1k,5k,10k,20k,100k]";
    if [ $# -lt 1 ]; then echo "$usage"; return; fi;

    local tmp=`mktemp -d`;
    local r=${2:-"1k,5k,10k,20k,100k"};
    r=`echo $r |  tr "," "\n" | str2num - | tr "\n" "," | sed "s/,$//"; `;
    hm mycat $1 > $tmp/a.pair
    cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 $c:${r%%,*} $tmp/a.pair $tmp/a.cool;
    cooler zoomify --balance -r$r $tmp/a.cool;
    cat $tmp/a.mcool
}


