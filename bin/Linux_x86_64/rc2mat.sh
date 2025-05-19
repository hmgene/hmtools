
rc2mat ()
{
    usage="$FUNCNAME <row,column,value> [min_per=0] [max_per=1]";
    if [ $# -lt 1 ]; then
        echo "$usage";
        return;
    fi;
    cat $1 | perl -e 'use strict; my ($m,$M)=('${2:-0}','${3:-1}');
        my %r=();
        my %c=();
        my $n=0;
        while(<STDIN>){chomp;my @d=split/\s+/,$_;
                $r{$d[0]}{$d[1]} += $d[2];
                $c{$d[1]}++;
                $n++;
        }
        my @cc=grep{ $c{$_}/$n >= $m && $c{$_}/$n <= $M } sort keys %c;
        print join("\t","rid",@cc),"\n";
        foreach my $i (sort keys %r){
                print join("\t",$i,map { defined $r{$i}{$_} ? $r{$i}{$_} : 0 } @cc),"\n";
        }
        '
}
rc2mat $@
