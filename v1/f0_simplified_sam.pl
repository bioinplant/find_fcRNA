open IN,"$ARGV[0]" or die "please input the sam file!";
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
open OUT, ">simplified_accepted_hits.sam";

while(<IN>){
	chomp;
	@aa=split/\t/,$_;
	@bb=split/ /,$aa[17];
	if($bb[0] eq "XF:Z:1"){#XF:Z:1 and XF:Z:2 provide the same information
		#push @out1, "$bb[1]\t$bb[2]\t$bb[3]";
		print OUT "$bb[1]\t$bb[2]\t$bb[3]\n";
		}
	}

##example of out
##Chr1-Chr1	11531	120m11513F30M
##Chr1-Chr1	11551	140m11513F10M
##Chr1-Chr2	12748	120m11352315F30M
##Chr1-Chr2	12727	96m11352315F54M
close IN;
close OUT;



