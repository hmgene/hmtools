
#!/bin/bash -l
BB=$BASEDIR/bin/`uname -sm | tr " " "_"`

str2num(){
    cat $1 | tr "[:lower:]" "[:upper:]" | numfmt --from=si
}

#LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-./}:$BB
ck(){
	echo $@ | tr " " "\n" | while read x;do
		[[ `command -v $x` ]] || { echo $x is not installed; }
	done 
}
groupby-sum(){
usage="$FUNCNAME <tsv> <key_columns> <value_column>"
if [ $# -lt 3 ];then echo "$usage";return;fi
	cat $1 | perl -e 'use strict; my @i=map{ $_ -1 } split /,/,"'$2'"; my $j='$3'-1; 
	my %r=(); while(<>){chomp;my@d=split/\s+/,$_;
		my $k=join("\t",map {$d[$_]} @i);
		$r{$k} += $d[$j];
	}
	foreach my $k (keys %r){
		print $k,"\t",$r{$k},"\n";
	}
	'
}
test-groupby-sum(){
echo "a	A	1
a	A	2
b	B	1
b	B	0.1" | groupby-sum - 1,2 3
}
groupby-max(){
usage="$FUNCNAME <tsv> <key_column> <value_column>"
if [ $# -lt 3 ];then echo "$usage";return;fi
	cat $1 | perl -e 'use strict; my $i='$2'-1; my $j='$3'-1; 
	my %r=(); while(<>){chomp;my@d=split/\s+/,$_;
		next if !defined $d[$i] || $d[$i] eq "";
		$r{$d[$i]}{$d[$j]} = $_;
	}
	foreach my $k (keys %r){
		my @x=sort {$b<=>$a} keys %{$r{$k}};
		print $r{$k}{$x[0]},"\n";
	}
	'
}
test-groupby-max(){
echo "a	1
a	2
b	1
b	0.1" | groupby-max - 1 2
}

check-required()
{
	ck bedtools samtools R
	
}

mycat(){
	if [ `file $1 | grep gzip | wc -l ` -gt 0 ];then
		gunzip -dc $1;
	else
		cat $1;
	fi
	
}
list-func(){
	ca declare -F | cut -d " " -f 3
}
jn(){
usage="$FUNCNAME <k> <c> <tsv> [<tsv>..]
join with summ operation 
";
if [ $# -lt 3 ];then echo "$usage";return;fi 

perl -e 'use strict; my %r=(); my $k='$1'-1;my $c='$2'-1;
	my $i=-1;
	foreach my $f (@ARGV){
		$i++;
		open( my $fh, "<", $f) or die "$!";
		while(<$fh>){chomp;my@d=split/\t/,$_;
			$r{$d[$k]}{$i}+=$d[$c];
		}
		close($fh);
	}
	foreach my $k (keys %r){
		print "$k";
		map { print "\t",defined $r{$k}{$_} ? $r{$k}{$_} : 0; } 0..$i;
		print "\n";
	}

' ${@:3}
}

jntest(){
	jn 1 2 <( echo -e "a\t1\nb\t2" ) <( echo -e "b\t11\nc\t22" )
}

repl()
{
    usage="$FUNCNAME <file> <key_value.txt> [replaced-only=0]
    multipe matches are delimitaated by |
    ";
    if [ $# -lt 2 ]; then echo "$usage"; return; fi;
    local tmp=$(mktemp -d );
    cat $2 > $tmp/a;
    cat $1 | perl -e 'use strict; my $only='${3:-0}';
                open(my $fh,"<","'$tmp/a'") or die $!;
		my %h=();map{chomp;my($k,$v)=split/\t/,$_; $h{$k}{$v}++;} <$fh>;
                close($fh);
                while(<STDIN>){chomp;  my @d=split/\t/,$_;
			my $hit=0;
			my @r=map {
				if( defined $h{$_} ){
					$hit++;
					$_=join("|",keys %{$h{$_}});
				}
				$_
			} @d; 
			if( $hit > 0 || $only == 0 ){
				print join("\t",@r),"\n";
			}
			
                }
        '
}
repl-test(){
echo "x	a	b
1	x	y
2	z	Z
3	x	z" | repl - <( echo -e "x\tX\ny\tY\ny\ty" ) 1
}


csv2xls(){
        cat $1| R --vanilla  -e 'tt=read.csv("stdin",header=T); write.table(file="'$2'",tt,quote=F,row.names=F,sep="\t")'
}


cutn(){
usage="$FUNCNAME <file> <column_name>[,<column_name>] [<head>=0]
        h=1 will include the head
"
if [ $# -lt 2 ];then echo "$usage";return;fi

local h=${3:-0};
cat $1 | perl -e 'use strict; my $h='${3:-0}'; my @cc=split/,/,"'$2'";
        my @c=();
        while(<STDIN>){chomp;$_=~s/\r//g; my@d=map{ $_ eq "" ? "NULL" : $_ } split/\t/,$_;
                if($#c < 0){
                        foreach my $j ( @cc){
                                map { if($j eq $d[$_]){ push @c,$_;} } 0..$#d;
                        }
                        if(!$h){ next;}
                }
                print join("\t",map {$d[$_]} @c),"\n";
        }
' ${@:1};
}


realpath(){
	echo `cd $( dirname $1 ); pwd -P`/$1
}

splitchrom(){
	cat $1  | awk -v odir=${2%\/} '/^>/{ gsub(">","",$1);o=odir"/"$1 ".fa"; print ">" $1 >> o; next}{ print >> o }' 
}

rose(){
        export PATH=$PATH:$BASEDIR/ROSE
        export PYTHONPATH=$BASEDIR/ROSE
	export WD=$BASEDIR/ROSE

        python2 $WD/ROSE_main.py $@
        #python2 $WD/ROSE_main.py -i $1 -r $2 -o $3 -g HG38  -t 2500
}
parallel(){
	$BB/$FUNCNAME $@
}
fastq-dump(){
	$BB/$FUNCNAME $@
}
fastp(){
	$BB/$FUNCNAME $@
}
fasta-get-markov(){
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BB
	$BB/$FUNCNAME $@
}

#samtools(){
#	$BB/$FUNCNAME $@
#}


crcmapper-rsem(){
        tail -n+2 $1 | cut -f 1,6 | perl -npe '$_=~s/ +/\t/g; $_=~s/\.\d+\t/\t/;'
}
#echo -e "#\nABC.3\t2\t3\t4\t5\t6" | crcmapper-rsem -

CRCMAPPER_HOME=/mnt/rstor/SOM_GENE_BEG33/users/hxk728/pr/CRCmapper
HG38_CHROM=/mnt/rstor/SOM_GENE_BEG33/users/hxk728/db/ucsc/hg38/chromosome/

crcmapper(){
	export WD=$BASEDIR/CRCmapper2
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BB
	python2 $WD/CRCmapper.py $@

#        local usage="$FUNCNAME <out_dir/name> <ROSE_output> <Peak.bed> <H3K27ac.bam> <chrom> [<TPM>]"
#        local odir=${1%\/*}
#        local nam=${1##*\/}
#        if [ $# -lt  5 ];then echo "$usage";return;fi
#        mkdir -p $odir
#        len=500; cof=33
#	local a=""
#        if [[ $# -gt 5 && -f $6 ]];then
#                a=" -a $6 ";
#        fi
#
#        python2 $WD/CRCmapper.py -o $odir/ -e $2 -s $3 -b $4 $a \
#        -g HG38 -f ${5%\/}/ -x $cof -l $len -n $nam
}
igvtools-bam2tdf(){
usage="$FUNCNAME <bam> <o.tdf> <ref.fa> [mem=8g]"
if [ $# -lt 3 ];then echo "$usage";return; fi
	java -Xmx${4:-8g} -jar $BASEDIR/IGVTools_2.3.98/igvtools.jar \
	count $1 $2 $3 
}


