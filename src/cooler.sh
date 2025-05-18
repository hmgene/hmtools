

mcool2bg(){
usage="$FUNCNAME <mcool> <resolution>";
if [ $# -lt 1 ];then echo "$usage";return;fi
local i=$1::resolutions/`echo ${2:-5k}| hm str2num`

python <( echo 'import cooler,sys
clr=cooler.Cooler("'$i'")
import src.neo
table,_=src.neo.cool2mar(clr)
table.to_csv(sys.stdout,sep="\t",index=False,header=False)
')
}

mcool2bp(){
usage="$FUNCNAME <mcool> <resolution>";
if [ $# -lt 2 ];then echo "$usage";return;fi 
    local r=`str2num ${2:-5k}`; 
    cooler dump --join  $1::resolutions/$r 
}

