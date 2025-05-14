



mcool2bg(){
usage="$FUNCNAME <mcool>";
if [ $# -lt 1 ];then echo "$usage";return;fi
python <( echo 'import cooler,sys
clr=cooler.Cooler("'$1'")
import src.neo
table,_=src.neo.cool2mar(clr)
table.to_csv(sys.stdout,sep="\t",index=False,header=False)
')
}

