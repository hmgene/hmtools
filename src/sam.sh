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


sam2bedpe () 
{ 
    usage="$FUNCNAME <sam> [<binsize=5k>] [min_span=1k]
        sam file must be paired -f0x1";
    if [ $# -lt 1 ]; then echo "$usage"; return; fi;
    local B=`echo ${2:-5k} | hm str2num -`;
    local S=`echo ${3:-1k} | hm str2num -`;

    cat $1 | perl -e ' my $min_span='$S'; my $B='$B';
    use strict;
    use warnings;

    sub cigar2len {
        my ($cigar) = @_; my $len = 0;
        while ($cigar =~ /(\d+)([MDN=X])/g) { $len += $1; }
        return $len;
    }
    while (<>) { chomp; next if /^@/;
        my @f = split(/\t/);

        my ($qname, $flag, $chr1, $pos1, $mapq, $cigar1, $rnext, $pnext, $tlen) = @f[0..8];
        my $chr2 = ($rnext eq "=") ? $chr1 : $rnext;
        next if $chr2 eq "*" || $pnext == 0;

        my $cigar_len1 = cigar2len($cigar1);
        my $start1 = $pos1 - 1;
        my $end1   = $start1 + $cigar_len1;
        my $strand1 = ($flag & 0x10) ? "-" : "+";
        my $strand2 = ($flag & 0x20) ? "-" : "+";

        # Find MC:Z: tag (mate CIGAR)
        my ($mc_tag) = grep { /^MC:Z:/ } @f;
        my $cigar2 = $mc_tag ? (split(/:/, $mc_tag, 3))[2] : "151M";  # fallback length
        my $cigar_len2 = cigar2len($cigar2);
        my $start2 = $pnext - 1;
        my $end2   = $start2 + $cigar_len2;

        # Consistent order to reduce redundancy
        if ($chr1 lt $chr2 || ($chr1 eq $chr2 && $start1 <= $start2)) {
            my $span = abs($end2 - $start1);
            print join("\t", $chr1, $start1, $end1, $chr2, $start2, $end2, $qname,$mapq, $strand1, $strand2), "\n";
        }
    }



'
}

sam2bed(){
usage="$FUNCNAME <sam|bam> [<binsz=100>]";
if [ $# -lt 1 ];then echo "$usage";return;fi 
    samtools view -F0x4 -bh $1 | bamToBed -split  
}
