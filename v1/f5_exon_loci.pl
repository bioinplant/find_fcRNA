use Getopt::Std;
use vars qw($opt_r);

# Usage
my $usage = "

f5_exon_loci.pl - exon loci of each sequences in temp4.txt .
      by  Qinjie Chu (Version: 1.0  Last Modified: 2016-12-02)

Usage: perl f5_exon_loci.pl [options]
 required:
  -r	refFlat file (simplfied).

";

getopts('r:');
die $usage unless $opt_r ;
$refflat = $opt_r if $opt_r;


open FLAT, "$refflat" or die "please input refFlat file!\n";
##format of refFlat file
#A1BG	NM_130786	chr19	-	58858171	58864865	58858387	58864803	8	58858171,58858718,58861735,58862756,58863648,58864293,58864657,58864769,	58858395,58859006,58862017,58863053,58863921,58864563,58864693,58864865,

while(<FLAT>){
	chomp;
	@cc=split/\t/,$_;
	$exons_num{$cc[0]}=$cc[8];
}
close FLAT;


open LEN, "fusion_len.txt";
open OUT, ">temp5.txt";

while(<LEN>){
	chomp;
	if(/>(.*)/){
		print OUT "$_\n";
		$fusion=$1;
		@aa=split/\|/,$fusion;
		$loc=0;
		
		@piars=();
		foreach $j(1..$exons_num{$aa[0]}){
			foreach $k(1..$exons_num{$aa[1]}){
				push @piars, "$aa[0]-$j|$aa[1]-$k";
				push @piars, "$aa[1]-$k|$aa[0]-$j";
				}
			if($j < $exons_num{$aa[0]}){
				$jj=$j+1;
				push @piars,"$aa[0]-$j|$jj";
				}
			}

		}
	else{
		@bb=split/\t/,$_;
		foreach $i(0..$#bb){
			$loc+=$bb[$i];
			print OUT "$loc\t$piars[$i]\n";
			}
		#print OUT "\n";
		}
	}
close LEN;
close OUT;


#format of temp5.txt
#>AAAS|CEP128
#258	AAAS-1|CEP128-1
#1363	CEP128-1|AAAS-1