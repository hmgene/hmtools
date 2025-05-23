
bedpe2bg-lr() 
{ 
    usage="$FUNCNAME <bedpe> <bin_size>";
    if [ $# -lt 2 ]; then echo "$usage"; return; fi;
    cat $1 | perl -e 'my %r=(); my $B2='$2'; my $B=int($B2/2);
	while(<STDIN>){chomp;my @d=split/\s+/,$_;
		next if ($#d < 5);
        my $v= $#d == 6 ? $d[$#d] : 1; ## handle 7 column bedpe and full bedpe 
        if ( $d[0] eq $d[3] ){
            $r{ join("\t",$d[0],int($d[1]/$B)) } += $v;
            $r{ join("\t",$d[3],int(($d[5]-1)/$B)) } -= $v;
        }
	}
	map{ my ($c,$s)=split/\t/,$_;
        print join("\t",$c,$s*$B,$s*$B + $B,$r{$_}),"\n";
    } keys %r;
    '
}
bedpe2bg-lr-test(){
    echo "c 1 2 c 4 5" | bedpe2bg-lr - 3
}

bedpe2bg() 
{ 
    usage="$FUNCNAME <bedpe> <bin_size>";
    if [ $# -lt 2 ]; then echo "$usage"; return; fi;
    cat $1 | perl -e 'my %r=(); my $B='$2';
	while(<STDIN>){chomp;my @d=split/\s+/,$_;
		next if ($#d < 5);
        my $v= $#d == 6 ? $d[$#d] : 1; ## handle 7 column bedpe and full bedpe 
        $r{ join("\t",$d[0],int($d[1]/$B)) } += $v;
        $r{ join("\t",$d[3],int(($d[5]-1)/$B)) } += $v;
	}
	map{ my ($c,$s)=split/\t/,$_;
        print join("\t",$c,$s*$B,$s*$B + $B,$r{$_}),"\n";
    } keys %r;
    '
}
bedpe2bg-test(){
    echo "c 1 2 c 4 5" | bedpe2bg - 3
}

bedpe-bin () 
{ 
    usage="$FUNCNAME <bedpe> <bin_size>";
    if [ $# -lt 2 ]; then echo "$usage"; return; fi;
    cat $1 | perl -e 'my %r=(); my $B='$2';
	while(<STDIN>){chomp;my @d=split/\s+/,$_;
		next if ($#d < 5);
        my $v= $#d == 6 ? $d[$#d] : 1; ## handle 7 column bedpe and full bedpe 
        $r{ join("\t",$d[0],int($d[1]/$B),$d[3],int(($d[5]-1)/$B)) } += $v;
	}
	map{ my ($c1,$s1,$c2,$e2)=split/\t/,$_;
        my ($s,$e)=($s1*$B,$e2*$B);
        print join("\t",$c1,$s,$s+$B,$c2,$e,$e+$B,$r{$_}),"\n";
    } keys %r;
    '
}

bedpe-bin-test(){
    echo "c 1 2 c 4 5" | bedpe-bin - 3
}

