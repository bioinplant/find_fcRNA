use Getopt::Std;
use vars qw( $opt_s $opt_l);

# Usage
my $usage = "

f1_fusion_site_from_simplified.pl - extract fusion sites from simplified sam file (tophat-fusion)
      by  Qinjie Chu (Version: 1.0  Last Modified: 2016-12-03)

Usage: perl f1_fusion_site_from_simplified.pl [options]
 required:
  -s	simplified sam file.
  -l	read length.

";

# command line processing.
getopts('s:l:');
die $usage unless $opt_s && $opt_l;

$sam = $opt_s if $opt_s;
$read_len = $opt_l if $opt_l;


open IN,"$sam" or die "please input the simplified sam file!";
##sam format
##1.reads name>>>>J00103:38:H3WYKBBXX:6:2102:4148:33123
##2.sam flag value>>>>129
##3.reference sequence name>>>>Chr1
##4.the left position of reference sequence(1-based)>>>>6909
##5.map quality>>>>50
##6.CIGAR>>>>114M
##7.reference name of the mate/next segment>>>>Chr10
##8.1-based leftmost position of the mate/next segment>>>>19040590
##9.observed Template LENgth>>>>0
##10.sequence>>>>CTACGTTTTGCCTTTGTTACCGATTCTACTAATTGTTGCTGATGAATGCGAATTCCGTGAAATGCCATCGCACAAGTGTTTTTATTCCTTTGTTTCCATTTCTACGGAGCAGGG
##11.sequence quality>>>>AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJFFJAAFFJFJJJJJJJJFJJJJJJJFJJJJJJJJJJJFJJFJJJJJ<JJJJJF7<AJ7<FJ)<
##12.unknown>>>>AS:i:0
##13.unknown>>>>XM:i:0
##14.unknown>>>>XO:i:0
##15.unknown>>>>XG:i:0
##16.unknown>>>>MD:Z:150
##17.unknown>>>>NM:i:0
##18.XF:Z:1 Chr1-Chr10 6909 114M19040738F36m CTACGTTTTGCCTTTGTTACCGATTCTACTAATTGTTGCTGATGAATGCGAATTCCGTGAAATGCCATCGCACAAGTGTTTTTATTCCTTTGTTTCCATTTCTACGGAGCAGGGGGGGGCAATAAGGTGGTGTTGGATGCCTGGTCGTTT AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJFFJAAFFJFJJJJJJJJFJJJJJJJFJJJJJJJJJJJFJJFJJJJJ<JJJJJF7<AJ7<FJ)<<))<-<-)-<<<7--<A7AF<A-A<A77A7<7FFF-
##19.XP:Z:Chr10 19040590 150M
##20.unknown>>>>NH:i:1

@out1=<IN>;

#while(<IN>){
#	chomp;
#	@aa=split/\t/,$_;
#	@bb=split/ /,$aa[17];
#	if($bb[0] eq "XF:Z:1"){#XF:Z:1 and XF:Z:2 provide the same information
#		push @out1, "$bb[1]\t$bb[2]\t$bb[3]";
#		#print OUT "$bb[1]\t$bb[2]\t$bb[3]\n";
#		}
#	}
#
##example of out1
##Chr1-Chr1	11531	120m11513F30M
##Chr1-Chr1	11551	140m11513F10M
##Chr1-Chr2	12748	120m11352315F30M
##Chr1-Chr2	12727	96m11352315F54M

foreach $out1(@out1){
	@cc=split/\t/,$out1;
	if($cc[2] =~ /(\d+)M(\d+)F(\d+)M/){
		$a=$cc[1]+$1-1;
		$b=$2+$3-1;
		push @out2, "$cc[0]\t$a\t$2\t$cc[2]";#only fusion sites is wanted
		}
	if($cc[2] =~ /(\d+)M(\d+)F(\d+)m/){
		$a=$cc[1]+$1-1;
		$b=$2-$3+1;
		push @out2, "$cc[0]\t$a\t$2\t$cc[2]";
		}
	if($cc[2] =~ /(\d+)m(\d+)F(\d+)M/){
		$a=$cc[1]-$1+1;
		$b=$2+$3-1;
		push @out2, "$cc[0]\t$a\t$2\t$cc[2]";
		}
	if($cc[2] =~ /(\d+)m(\d+)F(\d+)m/){
		$a=$cc[1]-$1+1;
		$b=$2-$3+1;
		push @out2, "$cc[0]\t$a\t$2\t$cc[2]";
		}
	}

#format of @out2
#CHROM	L-END	R-START	CIGAR

@out3=sort{(split/\t/,$a)[0] cmp (split/\t/,$b)[0] || (split/\t/,$a)[1] <=> (split/\t/,$b)[1] || (split/\t/,$a)[2] <=> (split/\t/,$b)[2]}@out2;


foreach (@out3){
	#print OUT "$_\n";
	@ee=split/\t/,$_;
	@ff=split/-/,$ee[0];
	
	##junction length over 10bp
	if($ee[3] =~ /(\d+)M(\d+)F(\d+)M/i){
		$len=$1;
		}
	
	##junction distance over 100kb between same chr; map length over 10bp
	if(($ff[0] eq $ff[1]) && ($ee[2]-$ee[1] > 100000) && ($len > 10 && $len+10 < $read_len)){
		push @out4, "$ee[0]\t$ee[1]\t$ee[2]";
		}
	
	##between different chrs, only map length over 10bp
	if(($ff[0] ne $ff[1]) && ($len > 10 && $len+10 < $read_len)){
		push @out4, "$ee[0]\t$ee[1]\t$ee[2]";
		}
	
	}

#example of @out4
#Chr1-Chr1	33738	23323804
#Chr1-Chr1	148105	23323653
#Chr1-Chr1	148105	23323653


open OUT, ">temp1.txt";
@out5=grep{++$count{$_}<2}@out4;
foreach (@out5){
	print OUT "$_\n";
	}
close OUT;

#example
#Chr1-Chr1	33738	23323804
#Chr1-Chr1	148105	23323653

