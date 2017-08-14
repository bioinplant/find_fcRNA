use Getopt::Std;
use vars qw($opt_b $opt_t $opt_l);

# Usage
my $usage = "

f6_from_bt1.pl - deal with the bowtie1 result.
      by  Qinjie Chu (Version: 1.0  Last Modified: 2017-01-05)

Usage: perl f6_from_bt1.pl [options]
 required:
  -b	bowtie1 result file.
  -t	temp5.txt that contains loci information of sequences in temp4.txt.
 optional:
  -l	the least length of bases that support one side of junctions. [ default = 10]

";

getopts('b:t:l:');
die $usage unless $opt_b && $opt_t;
$bowtieOut = $opt_b if $opt_b;
$temp = $opt_t if $opt_t;
$support_len = $opt_l ? $opt_l : 10 ;

open IN, "$bowtieOut" or die "please input bowtie1 out file.\n";
#format of bowtie out file
#SRR3239811.129 129 length=76	-	DHFR-DIAPH1	128093	CCAACATGATGATTGATGCAGCTAAGCTGCTTTCTGCTCTTTGTATTCTACCGCAGCCAGAGGACATGAATGAAAG	EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA	0	0:T>G,1:T>A
while(<IN>){
	chomp;
	@aa=split/\t/,$_;
	
	if($aa[7]=~/^(\d+):/ && $aa[7]!~/^(\d+):.+,(\d+):/){
				$map_start=$aa[3]+1;
				$map_end=$aa[3]+length($aa[4]);
				push @{$aa[2]}, $map_start;
				push @{$aa[2]}, $map_end;
				push @{$aa[2]}, $aa[0]; #read name
				push @{$aa[2]}, $aa[7]; #mismatch infomation
				push @{$aa[2]}, $aa[1]; #map strand
		}
	
	if($aa[7] =~ /^(\d+):.+,(\d+):/){
				$map_start=$aa[3]+1;
				$map_end=$aa[3]+length($aa[4]);
				push @{$aa[2]}, $map_start;
				push @{$aa[2]}, $map_end;
				push @{$aa[2]}, $aa[0];
				push @{$aa[2]}, $aa[7];
				push @{$aa[2]}, $aa[1];
		}
	
	if($aa[7] eq ""){
				$map_start=$aa[3]+1;
				$map_end=$aa[3]+length($aa[4]);
				push @{$aa[2]}, $map_start;
				push @{$aa[2]}, $map_end;
				push @{$aa[2]}, $aa[0];
				push @{$aa[2]}, " ";
				push @{$aa[2]}, $aa[1];
		}
	
	}
close IN;

open TEMP, "$temp" or die "please input temp5.txt that contains loci information of sequences in temp4.txt .\n";
#format of temp6_exons.txt
#>AAAS|CEP128
#258	AAAS-1|CEP128-1
#1363	CEP128-1|AAAS-1

open OUT, ">temp6.txt";
while(<TEMP>){
	chomp;
	if(/>(.*)/){
		print OUT "$_\n";
		$fusion=$1;
		}
	else{
		@cc=split/\t/,$_;

			for($j=0;$j<=$#{"$fusion"};$j=$j+5){
				if(${"$fusion"}[$j+4] eq '+'){
					#${"$fusion"}[$j] is start position, 1-based
					#${"$fusion"}[$j+1] is end position, 1-based
					#${"$fusion"}[$j+2] is read name
					#${"$fusion"}[$j+3] is read mismatch information
					#${"$fusion"}[$j+4] is map strand
					#接口左右3-10bp不能有错配
					#接口左右0-2bp可以有错配
					if(${"$fusion"}[$j+3]=~/^(\d+):/ && (${"$fusion"}[$j+3])!~/^(\d+):.+,(\d+):/){ #if there is only one mismatch
						$loci=$1;
						if((${"$fusion"}[$j]+$loci+9<$cc[0] || ${"$fusion"}[$j]+$loci-10>$cc[0] || ( ${"$fusion"}[$j]+$loci <= $cc[0] && ${"$fusion"}[$j]+$loci+2 >= $cc[0]) || ( ${"$fusion"}[$j]+$loci >= $cc[0] && ${"$fusion"}[$j]+$loci-3 <= $cc[0])) && ${"$fusion"}[$j] - 1 + $support_len <$cc[0] && (${"$fusion"}[$j+1]-$support_len >$cc[0]) && ($cc[1] !~ /\d+\|\d+/)){
							print OUT "${\"$fusion\"}[$j]\t${\"$fusion\"}[$j+1]\t$cc[1]\t${\"$fusion\"}[$j+2]\n";
							}
						}
				
					if(${"$fusion"}[$j+3] =~ /^(\d+):.+,(\d+):/){ #if there is two mismatches
						$loci1=$1;
						$loci2=$2;
						if((${"$fusion"}[$j]+$loci1+9<$cc[0] || ${"$fusion"}[$j]+$loci1-10>$cc[0]  || ( ${"$fusion"}[$j]+$loci1 <= $cc[0] && ${"$fusion"}[$j]+$loci1+2 >= $cc[0]) || ( ${"$fusion"}[$j]+$loci1 >= $cc[0] && ${"$fusion"}[$j]+$loci1-3 <= $cc[0])) && (${"$fusion"}[$j]+$loci2+9<$cc[0] || ${"$fusion"}[$j]+$loci2-10>$cc[0] || ( ${"$fusion"}[$j]+$loci2 <= $cc[0] && ${"$fusion"}[$j]+$loci2+2 >= $cc[0]) || ( ${"$fusion"}[$j]+$loci2 >= $cc[0] && ${"$fusion"}[$j]+$loci2-3 <= $cc[0])) && ${"$fusion"}[$j] - 1 +$support_len<$cc[0] && (${"$fusion"}[$j+1]-$support_len>$cc[0]) && ($cc[1] !~ /\d+\|\d+/)){
							print OUT "${\"$fusion\"}[$j]\t${\"$fusion\"}[$j+1]\t$cc[1]\t${\"$fusion\"}[$j+2]\n";
							}
						}
				
					if(${"$fusion"}[$j+3] eq " "){ #if there is no mismatch
						if(${"$fusion"}[$j] - 1 +$support_len <$cc[0] && (${"$fusion"}[$j+1]- $support_len>$cc[0]) && ($cc[1] !~ /\d+\|\d+/)){
							print OUT "${\"$fusion\"}[$j]\t${\"$fusion\"}[$j+1]\t$cc[1]\t${\"$fusion\"}[$j+2]\n";
							}
						}
				}

				if(${"$fusion"}[$j+4] eq '-'){
					$read_len = ${"$fusion"}[$j+1] - ${"$fusion"}[$j] + 1;
					if(${"$fusion"}[$j+3]=~/^(\d+):/ && (${"$fusion"}[$j+3])!~/^(\d+):.+,(\d+):/){ #if there is only one mismatch
						$loci = $read_len - $1 - 1;
						if((${"$fusion"}[$j]+$loci+9<$cc[0] || ${"$fusion"}[$j]+$loci-10>$cc[0] || ( ${"$fusion"}[$j]+$loci <= $cc[0] && ${"$fusion"}[$j]+$loci+2 >= $cc[0]) || ( ${"$fusion"}[$j]+$loci >= $cc[0] && ${"$fusion"}[$j]+$loci-3 <= $cc[0])) && ${"$fusion"}[$j]-1+$support_len<$cc[0] && (${"$fusion"}[$j+1]-$support_len>$cc[0]) && ($cc[1] !~ /\d+\|\d+/)){
							print OUT "${\"$fusion\"}[$j]\t${\"$fusion\"}[$j+1]\t$cc[1]\t${\"$fusion\"}[$j+2]\n";
							}
						}
				
					if(${"$fusion"}[$j+3] =~ /^(\d+):.+,(\d+):/){ #if there is two mismatches
						$loci1 = $read_len - $1 - 1;
						$loci2 = $read_len - $2 - 1;
						if((${"$fusion"}[$j]+$loci1+9<$cc[0] || ${"$fusion"}[$j]+$loci1-10>$cc[0] || ( ${"$fusion"}[$j]+$loci1 <= $cc[0] && ${"$fusion"}[$j]+$loci1+2 >= $cc[0]) || ( ${"$fusion"}[$j]+$loci1 >= $cc[0] && ${"$fusion"}[$j]+$loci1-3 <= $cc[0])) && (${"$fusion"}[$j]+$loci2+9<$cc[0] || ${"$fusion"}[$j]+$loci2-10>$cc[0] || ( ${"$fusion"}[$j]+$loci2 <= $cc[0] && ${"$fusion"}[$j]+$loci2+2 >= $cc[0]) || ( ${"$fusion"}[$j]+$loci2 >= $cc[0] && ${"$fusion"}[$j]+$loci2-3 <= $cc[0])) && ${"$fusion"}[$j]-1+$support_len<$cc[0] && (${"$fusion"}[$j+1]-$support_len>$cc[0]) && ($cc[1] !~ /\d+\|\d+/)){
							print OUT "${\"$fusion\"}[$j]\t${\"$fusion\"}[$j+1]\t$cc[1]\t${\"$fusion\"}[$j+2]\n";
							}
						}
				
					if(${"$fusion"}[$j+3] eq " "){ #if there is no mismatch
						if(${"$fusion"}[$j]-1+$support_len<$cc[0] && (${"$fusion"}[$j+1]-$support_len>$cc[0]) && ($cc[1] !~ /\d+\|\d+/)){
							print OUT "${\"$fusion\"}[$j]\t${\"$fusion\"}[$j+1]\t$cc[1]\t${\"$fusion\"}[$j+2]\n";
							}
						}
					
					}
				
				}

		}
	}
close TEMP;
close OUT;