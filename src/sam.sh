
sam2bed(){
usage="$FUNCNAME <sam|bam> [<binsz=100>]";
if [ $# -lt 1 ];then echo "$usage";return;fi 
    samtools view -F0x4 -bh $1 | bamToBed -split  
}
