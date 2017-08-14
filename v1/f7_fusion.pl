use warnings;
use Getopt::Std;
use vars qw($opt_r $opt_t $opt_s);


# Usage
my $usage = "

f7_fusion.pl - final step to identify fusion circRNAs.
      by  Qinjie Chu (Version: 1.0  Last Modified: 2017-01-05)

Usage: perl f7_fusion.pl [options]
 required:
  -r	refFlat file.
  -t	temp6.txt that contains bowtie1 map information.
 optional:
  -s	the least number of support reads at each junctions. [ default = 2 ]

";

getopts('r:t:s:');
die $usage unless $opt_r && $opt_t;
$refflat = $opt_r if $opt_r;
$temp = $opt_t if $opt_t;
$supportreads = $opt_s ? $opt_s : 2;

open FLAT, "$refflat" or die "please input refFlat file!\n";
##format of refFlat file
#A1BG	NM_130786	chr19	-	58858171	58864865	58858387	58864803	8	58858171,58858718,58861735,58862756,58863648,58864293,58864657,58864769,	58858395,58859006,58862017,58863053,58863921,58864563,58864693,58864865,
while(<FLAT>){
	chomp;
	@ff=split/\t/,$_;
	$ori{$ff[0]}=$ff[3];
}
close FLAT;

open IN, "$temp" or die "please input temp6.txt !\n";
#format of temp6.txt
#>C15F1.8|Y55F3AM.21
#5618	5718	Y55F3AM.21-7|C15F1.8-2 reads_name
open OUT, ">temp7.txt";
open OUT2, ">temp7_name.txt";

@ab=();
@ba=();
%count=();
%count2=();

while(<IN>){
	chomp;
	if(/>(.+)/){
		
		@aa=split/\|/,$1;
		
		@ab=grep {++$count{$_} >= $supportreads;}@ab;
		@ba=grep {++$count{$_} >= $supportreads;}@ba;
		@unique_ab=grep {++$count2{$_} < 2;}@ab;
		@unique_ba=grep {++$count2{$_} < 2;}@ba;
		
		
		foreach $i(@unique_ab){
			#$i format: geneA	geneB	fusion_exon_of_geneA	fusion_exon_of_geneB
			@cc=split/\t/,$i;
			
			foreach $j(@unique_ba){
				#$j format: geneB	geneA	fusion_exon_of_geneB	fusion_exon_of_geneA
				@dd=split/\t/,$j;
				
				if(($ori{$cc[0]} eq "+") && ($ori{$cc[1]} eq "+") && (($dd[2]>=$cc[3] && ($dd[3]<=$cc[2])) || ($cc[2]>=$dd[3] && ($cc[3]<=$dd[2])))){
					foreach (@{"$cc[0]\-$cc[2]|$cc[1]\-$cc[3]"}){$readname1.="$_\,";}
					foreach (@{"$dd[0]\-$dd[2]|$dd[1]\-$dd[3]"}){$readname2.="$_\,";}
					push @final, "$cc[0]|$cc[1]\t$cc[2]\t$cc[3]\t$count{$i}reads\t$readname1\n$dd[0]|$dd[1]\t$dd[2]\t$dd[3]\t$count{$j}reads\t$readname2\n";
					push @name, "$cc[0]|$cc[1]";
					$readname1="";
					$readname2="";
					}
				
				if(($ori{$cc[0]} eq "+") && ($ori{$cc[1]} eq "-") && (($dd[2]<=$cc[3] && ($dd[3]<=$cc[2])) || ($cc[2]>=$dd[3] && ($cc[3]>=$dd[2])))){
					foreach (@{"$cc[0]\-$cc[2]|$cc[1]\-$cc[3]"}){$readname1.="$_\,";}
					foreach (@{"$dd[0]\-$dd[2]|$dd[1]\-$dd[3]"}){$readname2.="$_\,";}
					push @final, "$cc[0]|$cc[1]\t$cc[2]\t$cc[3]\t$count{$i}reads\t$readname1\n$dd[0]|$dd[1]\t$dd[2]\t$dd[3]\t$count{$j}reads\t$readname2\n";
					push @name, "$cc[0]|$cc[1]";
					$readname1="";
					$readname2="";
					}
				
				if(($ori{$cc[0]} eq "-") && ($ori{$cc[1]} eq "+") && (($dd[2]>=$cc[3] && ($dd[3]>=$cc[2])) || ($cc[2]<=$dd[3] && ($cc[3]<=$dd[2])))){
					foreach (@{"$cc[0]\-$cc[2]|$cc[1]\-$cc[3]"}){$readname1.="$_\,";}
					foreach (@{"$dd[0]\-$dd[2]|$dd[1]\-$dd[3]"}){$readname2.="$_\,";}
					push @final, "$cc[0]|$cc[1]\t$cc[2]\t$cc[3]\t$count{$i}reads\t$readname1\n$dd[0]|$dd[1]\t$dd[2]\t$dd[3]\t$count{$j}reads\t$readname2\n";
					push @name, "$cc[0]|$cc[1]";
					$readname1="";
					$readname2="";
					}
				
				if(($ori{$cc[0]} eq "-") && ($ori{$cc[1]} eq "-") && (($dd[2]<=$cc[3] && ($dd[3]>=$cc[2])) || ($cc[2]<=$dd[3] && ($cc[3]>=$dd[2])))){
					foreach (@{"$cc[0]\-$cc[2]|$cc[1]\-$cc[3]"}){$readname1.="$_\,";}
					foreach (@{"$dd[0]\-$dd[2]|$dd[1]\-$dd[3]"}){$readname2.="$_\,";}
					push @final, "$cc[0]|$cc[1]\t$cc[2]\t$cc[3]\t$count{$i}reads\t$readname1\n$dd[0]|$dd[1]\t$dd[2]\t$dd[3]\t$count{$j}reads\t$readname2\n";
					push @name, "$cc[0]|$cc[1]";
					$readname1="";
					$readname2="";
					}
				}
			}

		@unique_ba=();
		@unique_ab=();
		@ab=();
		@ba=();
		%count=();
		%count2=();
		}
	else{
		
		@bb=split/\t/,$_;
		if($bb[2] =~ /$aa[0]-(\d+)\|$aa[1]-(\d+)/){
			push @ab, "$aa[0]\t$aa[1]\t$1\t$2";
			#geneA	geneB	fusion_exon_of_geneA	fusion_exon_of_geneB
			push @{$bb[2]}, $bb[3];
			}
		if($bb[2] =~ /$aa[1]-(\d+)\|$aa[0]-(\d+)/){
			push @ba, "$aa[1]\t$aa[0]\t$1\t$2";
			#geneB	geneA	fusion_exon_of_geneB	fusion_exon_of_geneA
			push @{$bb[2]}, $bb[3];
			}

		}
	}

		@ab=grep {++$count{$_} >= $supportreads;}@ab;
		@ba=grep {++$count{$_} >= $supportreads;}@ba;
		@unique_ab=grep {++$count2{$_} < 2;}@ab;
		@unique_ba=grep {++$count2{$_} < 2;}@ba;

foreach $i(@unique_ab){
			@cc=split/\t/,$i;
			
			foreach $j(@unique_ba){
				@dd=split/\t/,$j;
				
				if(($ori{$cc[0]} eq "+") && ($ori{$cc[1]} eq "+") && (($dd[2]>=$cc[3] && ($dd[3]<=$cc[2])) || ($cc[2]>=$dd[3] && ($cc[3]<=$dd[2])))){
					foreach (@{"$cc[0]\-$cc[2]|$cc[1]\-$cc[3]"}){$readname1.="$_\,";}
					foreach (@{"$dd[0]\-$dd[2]|$dd[1]\-$dd[3]"}){$readname2.="$_\,";}
					push @final, "$cc[0]|$cc[1]\t$cc[2]\t$cc[3]\t$count{$i}reads\t$readname1\n$dd[0]|$dd[1]\t$dd[2]\t$dd[3]\t$count{$j}reads\t$readname2\n";
					push @name, "$cc[0]|$cc[1]";
					$readname1="";
					$readname2="";
					}
				
				if(($ori{$cc[0]} eq "+") && ($ori{$cc[1]} eq "-") && (($dd[2]<=$cc[3] && ($dd[3]<=$cc[2])) || ($cc[2]>=$dd[3] && ($cc[3]>=$dd[2])))){
					foreach (@{"$cc[0]\-$cc[2]|$cc[1]\-$cc[3]"}){$readname1.="$_\,";}
					foreach (@{"$dd[0]\-$dd[2]|$dd[1]\-$dd[3]"}){$readname2.="$_\,";}
					push @final, "$cc[0]|$cc[1]\t$cc[2]\t$cc[3]\t$count{$i}reads\t$readname1\n$dd[0]|$dd[1]\t$dd[2]\t$dd[3]\t$count{$j}reads\t$readname2\n";
					push @name, "$cc[0]|$cc[1]";
					$readname1="";
					$readname2="";
					}
				
				if(($ori{$cc[0]} eq "-") && ($ori{$cc[1]} eq "+") && (($dd[2]>=$cc[3] && ($dd[3]>=$cc[2])) || ($cc[2]<=$dd[3] && ($cc[3]<=$dd[2])))){
					foreach (@{"$cc[0]\-$cc[2]|$cc[1]\-$cc[3]"}){$readname1.="$_\,";}
					foreach (@{"$dd[0]\-$dd[2]|$dd[1]\-$dd[3]"}){$readname2.="$_\,";}
					push @final, "$cc[0]|$cc[1]\t$cc[2]\t$cc[3]\t$count{$i}reads\t$readname1\n$dd[0]|$dd[1]\t$dd[2]\t$dd[3]\t$count{$j}reads\t$readname2\n";
					push @name, "$cc[0]|$cc[1]";
					$readname1="";
					$readname2="";
					}
				
				if(($ori{$cc[0]} eq "-") && ($ori{$cc[1]} eq "-") && (($dd[2]<=$cc[3] && ($dd[3]>=$cc[2])) || ($cc[2]<=$dd[3] && ($cc[3]>=$dd[2])))){
					foreach (@{"$cc[0]\-$cc[2]|$cc[1]\-$cc[3]"}){$readname1.="$_\,";}
					foreach (@{"$dd[0]\-$dd[2]|$dd[1]\-$dd[3]"}){$readname2.="$_\,";}
					push @final, "$cc[0]|$cc[1]\t$cc[2]\t$cc[3]\t$count{$i}reads\t$readname1\n$dd[0]|$dd[1]\t$dd[2]\t$dd[3]\t$count{$j}reads\t$readname2\n";
					push @name, "$cc[0]|$cc[1]";
					$readname1="";
					$readname2="";
					}
				}
			}
	
	
	$dup{$_}++ for @final;
	foreach $key(sort keys %dup){
		print OUT "$key\n";
		}
	
	$dup2{$_}++ for @name ;
	foreach $key(sort keys %dup2){
		print OUT2 "$key\n";
		}