#!/bin/bash
#SBATCH -c 32
#SBATCH --mem=96g


hg38=/mnt/vstor/SOM_GENE_BEG33/ref_files/hg38/bwa_indices/GRChg38.d1.vd1.fa
ecoli=/mnt/vstor/SOM_GENE_BEG33/ref_files/Escherichia_coli_K_12_MG1655/Sequence/BWAIndex/genome.fa
dm3=/mnt/vstor/SOM_GENE_BEG33/ref_files/dm3/dm3.fa
mm10=/mnt/vstor/SOM_GENE_BEG33/ref_files/mm10/bwa_indices/mm10.fa
ercc=/mnt/vstor/SOM_GENE_BEG33/RNA_seq/hg19/ref/RSEM_ERCC/ERCC92.fa
genomes=( hg38 ecoli dm3 mm10 ercc)
export ${genomes[@]};


#sample=SaOS-LM7_invivo_w8_r1_ATAC_120723
input=( 
/mnt/vstor/SOM_GENE_BEG33/data/120723/SaOS-LM7_invivo_w8_r1_ATAC_120723_R1.fastq.gz  
/mnt/vstor/SOM_GENE_BEG33/data/120723/SaOS-LM7_invivo_w8_r3_ATAC_120723_R1.fastq.gz
)
output=o
p1(){
	s=${1##*/};s=${s%\_R1.fastq.gz};
	o=$2/$s/fq
	mkdir -p $o
	if [ -s $1 -a ! -s $o/R1.aa ];then
		zcat $1 | split -l40000000 - $o/R1.
	fi
	if [ -s ${1/_R1/_R2} -a ! -s $o/R2.aa ];then
		zcat ${1/_R1/_R2} | split -l40000000 - $o/R2.
	fi
}
p2(){
	i1=$1;i2=${i1/\/R1\./\/R2\.};r=$2;o=${i1/fq\/R1./bwa\/}@$r.bam
	if [ -s $i1 -a -s $i2 -a ${!r} -a ! -s $o ];then
		mkdir -p ${o%\/*}
		bwa mem -M -t 8 ${!r} $i1 $i2 | samtools view -bhq 30  | samtools sort -T $o.tmp - > $o
		bam-score $o | gzip -c > $o.id.gz
	fi
}

bam-score(){
	samtools view -F0x4 $1 | perl -e 'use strict; my %h=();
        while(<STDIN>){chomp; my @d=split/\t/,$_;
		$h{$d[0]}{$d[2]}{$d[1] & 0x40 ? 1 : 2 }{ join(",",@d[3..5]) } ++; ## pos, qc, cigar
        }
	sub sc{ my ($x)=@_; my $r=0;while($x=~/(\d+)M/g){$r+=$1;}; return $r; }
	sub score{ my ($x,$y)=@_;
		my $type= scalar @$x > 0 && scalar @$y ? 2 : 1;
		my $dist=0;
		my $score=0;
		if($type == 2){
			foreach my $i (@$x){
				my ($ix,$is,$ic)=split/,/,$i;
				foreach my $j (@$y){
					my ($jx,$js,$jc)=split/,/,$j;
					if( $score < sc($ic.$jc) ){
						$score = sc($ic.$jc);
						$dist = abs($ix-$jx);
					}
				}
			}
		}else{
			foreach my $i (@$x,@$y){
				my ($ix,$is,$ic)=split/,/,$i;
				if( $score < sc($ic) ){
					$score = sc($ic);
					$dist = 0;
				}
			}
		}
		return $type.$score;
	}

	foreach my $k (keys %h){
		my $hscore=0;
		foreach my $c (keys %{$h{$k}}){
			my @x=keys %{$h{$k}{$c}{1}};
			my @y=keys %{$h{$k}{$c}{2}};
			my $s=score(\@x,\@y); 
			$hscore = $hscore > $s ? $hscore : $s;
		}
		print "$k\t$hscore\n";

	}
' 
}

p3(){
	s=${1#*\@};s=${s%.bam*}; ## this species
	v=( `ls ${1%\@*}*.bam.id.gz | grep -v $s` ); ## other species
	o=${1/bwa/bwa_flt}
	mkdir -p ${o%\/*}
	echo "$1 => $o";
	gunzip -dc "${v[@]}" | perl -e 'use strict; my %r=(); 
	map {chomp;my($x,$y)=split/\t/,$_; $r{$x}=defined $r{$x} && $r{$x} < $y ? $r{$x} : $y;  } <STDIN>;
	map{ chomp;my($x,$y)=split/\t/,$_;
		if( ! defined $r{$x} || $r{$x} < $y ){ 
			print $_,"\n"; 
		}
	} `gunzip -dc '"${1/\@/\\@}"'`;
	'  | gzip -c > $o
	i=${1%.id.gz};
	i2=${1/bwa/bwa_flt};
	o=${i2%.id.gz}
	echo "$i => $o"
	{ 
		samtools view -H $i
		samtools view $i | perl -e 'use strict; my %r=();
		map {chomp; $r{$_}++;} `gunzip -dc '${i2/\@/\\@}' | cut -f 1`;
		while(<>){chomp;my@d=split/\t/,$_;
			if(defined $r{$d[0]} ){ print $_,"\n";}
		}
		'  
	} | samtools view -bh > $o 


#samtools view -hF0x4 $2 | perl -e 'use strict; my %h=();
#        my %h=map {chomp;my @d=split/\t/,$_; $d[0]=>1 } `cat '$4.ids'`;
#        while(<STDIN>){chomp; my @d=split/\t/,$_;
#                if($_=~/^@/ || defined $h{$d[0]} ){ print $_,"\n";}
#        }
#' | samtools view -hb > $4
#samtools index $4; samtools flagstat $4 > ${4%.bam}.flagstat.txt
}
p4(){
	x=( `ls $1/*\@$2.bam` );
	n=${x[0]%\/bwa_flt*};n=${n##*/};
	o=${x[0]%\/bwa_flt*}/$n\@$2.bam
	echo "${x[0]}...=>$o";
	{
		samtools view -H ${x[0]};
		for f in ${x[@]};do
			samtools view $f
		done
	} | samtools view -bh | samtools sort  -T $o.tmp > $o
	samtools index $o
	samtools flagstat $o > $o.flagstat.txt
}
export -f p1 p2 p3 p4 bam-score 
#parallel -j1 time p1 {} $output ::: ${input[0]};
#parallel -j4 p2 {} ::: $output/*/*/R1.* ::: ${genomes[@]};
parallel -j8 p3 {} ::: $output/*/bwa/*.bam.id.gz  
parallel -j8 p4 {} ::: $output/*/bwa_flt ::: ${genomes[@]};

