bedw(){
usage="$FUNCNAME <bed> <window>";
if [ $# -lt 2 ];then echo "$usage";return;fi;
    w=`echo $2 | hm str2num -`
    cat | awk -v OFS="\t" -v w=$w '{ $2=$2-w;$3=$3+w;print $0;}'
}

bed5p(){
if [ $# -lt 1 ];then echo "$usage";return;fi;
    cat $1 | perl -ne 'chomp;my($c,$s,$e,$n,$v,$t)=split/\t/,$_;my $p= $t eq "+" ? $s : $e -1;  print join("\t",$c,$p,$p+1,$n,$v,$t),"\n";'
}


bg2bw () 
{ 
    usage="$FUNCNAME <bg> [<genome=hg38>]";
    local tmpd=$( mktemp -d );
    hm ucsc-chromsize ${2:-hg38} > $tmpd/b;
    cat $1 | perl -e 'use strict;
		my %r=();
		map {chomp; my($k,$v)=split/\t/,$_; $r{$k}=$v; } `cat '$tmpd/b'`;
		while(<>){chomp; my @d=split/\t/,$_;
			$d[1] = $d[1] < 0 ? 0 : $d[1];
			if( $d[2]  < $r{$d[0]} ){	
				print join("\t",@d),"\n";
			}
		}
	' | sort -k1,1 -k2,3n > $tmpd/a;
    `ca home`/bin/bedGraphToBigWig $tmpd/a $tmpd/b $tmpd/o;
    cat $tmpd/o
}


bw2bg () 
{ 
    usage="$FUNCNAME <bw> [options]";
    if [ $# -lt 1 ]; then
        echo "$usage";
        return;
    fi;
    local tmpd=$( mktemp -d );
    hm mycat $1 > $tmpd/a;
    $HMHOME/bin/bigWigToBedGraph $tmpd/a $tmpd/o ${@:2};
    cat $tmpd/o
}

bed2bg () 
{ 
    usage="$FUNCNAME <bed> [<binsize>=100]
";
    if [ $# -lt 1 ]; then
        echo "$usage";
        return;
    fi;
    cat $1 | perl -e 'use strict; my $B='${2:-100}';
    my %r=(); 
    while(<>){ chomp; my @d=split/\s+/,$_;
        my ($s,$e)=(int($d[1]/$B),int(($d[2]-1)/$B));
        for(my$i=$s;$i<=$e;$i++){
            $r{$d[0]}{$i} ++;
        }
    }
    foreach my $c (keys %r){
        map { print join("\t",$c,$_*$B,($_+1)*$B,$r{$c}{$_}),"\n"; } sort {$a<=>$b} keys %{$r{$c}};
    }
        
'
}


bed2bg-test(){
echo "--input--"
echo "c 0 11
c 10 20" 
echo "--output--"
echo "c 0 11
c 10 20" | bed2bg - 10
}

