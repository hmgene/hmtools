
ucsc-refflat2bed12(){
cat $1 | perl -F'\t' -lane '
  $name = $F[1];
  $chrom = $F[2];
  $strand = $F[3];
  $txStart = $F[4];
  $txEnd = $F[5];
  $cdsStart = $F[6];
  $cdsEnd = $F[7];
  $exonCount = $F[8];
  @starts = split(",", $F[9]);
  @ends = split(",", $F[10]);
  @blockSizes = map { $ends[$_] - $starts[$_] } (0..$#starts);
  @blockStarts = map { $starts[$_] - $txStart } (0..$#starts);
  print join("\t", $chrom, $txStart, $txEnd, $name, 0, $strand, $cdsStart, $cdsEnd, 0,
              $exonCount, join(",", @blockSizes), join(",", @blockStarts));
' 
}

ucsc-dn(){
usage="$FUNCNAME <genome> <file>"; if [ $# -lt 2 ];then echo "$usage";return;fi
    echo "downloading .."
    wget -O $2 https://hgdownload.soe.ucsc.edu/goldenPath/$1/database/$2
}
ucsc-refflat(){
usage="$FUNCNAME <genome=hg38>"; if [ $# -lt 1 ];then echo "$usage";return;fi
    f=$HMHOME/data/ucsc/$1/refFlat.txt.gz
    if [ -f $f ];then 
        hm mycat $f
    else
        hm ucsc-dn $1 refFlat.txt.gz
        mkdir -p ${f%/*}; echo "mv refFlat.txt.gz ${f%/*}"
        mv refFlat.txt.gz ${f%/*}
        hm mycat $f
    fi
}
ucsc-chromsize(){
usage="$FUNCNAME <genome=hg38>"; if [ $# -lt 1 ];then echo "$usage";return;fi
    f=$HMHOME/data/ucsc/$1/chromInfo.txt.gz
    if [ -f $f ];then 
        hm mycat $f
    else
        hm ucsc-dn $1 chromInfo.txt.gz
        mkdir -p ${f%/*}; echo "mv chromInfo.txt.gz ${f%/*}"
        mv chromInfo.txt.gz ${f%/*}
        hm mycat $f 
    fi
}

