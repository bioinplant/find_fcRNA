use Getopt::Std;
use vars qw( $opt_g $opt_r $opt_t );

# Usage
my $usage = "

f4_transcript_exon_seq_rc.pl - extract exon sequences of transcript pairs.
      by  Qinjie Chu (Version: 1.0  Last Modified: 2016-12-29)

Usage: perl f4_transcript_exon_seq_rc.pl [options]
 required:
  -g	genome file.
  -r	refFlat file (simplfied).
  -t	temp3.txt that contains the fusion transcript pairs.

";

getopts('g:r:t:');
die $usage unless $opt_g && $opt_r && $opt_t ;

$genome = $opt_g if $opt_g;
$refflat = $opt_r if $opt_r;
$temp = $opt_t if $opt_t;


open GENOME, "$genome" or die "please input genome file!\n";
while(<GENOME>){
	chomp;
	if(/>(.*)/){
		$name=$1;
		}
	else{
		$ha{$name}.=$_;
		}
	}
close GENOME;


open FLAT, "$refflat" or die "please input refFlat file!\n";
##format of refFlat file
#A1BG	NM_130786	chr19	-	58858171	58864865	58858387	58864803	8	58858171,58858718,58861735,58862756,58863648,58864293,58864657,58864769,	58858395,58859006,58862017,58863053,58863921,58864563,58864693,58864865,

while(<FLAT>){
	chomp;
	@aa=split/\t/,$_;
	
	@start1=split/,/,$aa[9];
	@end1=split/,/,$aa[10];
	
			foreach $i(0..$#start1){
				$exons=$i+1;
				##start position in refFlat.txt downloaded from UCSC is 0-based
				##end position in refFlat.txt downloaded from UCSC is 1-based
				if(exists $ha{$aa[2]}){
					#$out= uc(substr($ha{$aa[2]},$start1[$i],($end1[$i]-$start1[$i]+1)));
					$out= uc(substr($ha{$aa[2]},$start1[$i],($end1[$i]-$start1[$i])));
					if($aa[3] eq '+'){
						$seq=$out;
						push @{$aa[0]}, "$aa[0]-$exons\t$seq";
						}
					if($aa[3] eq '-'){
						$seq=reverse($out);
						$seq=~tr/ACGTacgt/TGCAtgca/;
						push @{$aa[0]}, "$aa[0]-$exons\t$seq";
						}
					}
			}
	
	}
close FLAT;

open IN, "$temp" or die "please input temp3.txt that contains the fusion transcript pairs.\n";
##format of temp3.txt
#CAMTA1	MTF1
#ENO1	EDARADD
#ENO1	EDARADD
open OUT, ">temp4.txt";
while(<IN>){
	chomp;
	if(!exists $hapair{$_}){
		$hapair{$_}=1;
		@bb=split/\t/,$_;
		print OUT ">$bb[0]|$bb[1]\n";
		push @fusion, "$bb[0]|$bb[1]";
		foreach $left(@{$bb[0]}){
			@cc=split/\t/,$left;
			$leftlen=length($cc[1]);
			print OUT "$cc[1]";
			push @{"$bb[0]|$bb[1]"}, "$leftlen";
				foreach $right(@{$bb[1]}){
					@dd=split/\t/,$right;
					$rightlen=length($dd[1]);
					print OUT "$dd[1]$cc[1]";
					push @{"$bb[0]|$bb[1]"}, "$rightlen";
					push @{"$bb[0]|$bb[1]"}, "$leftlen";
					}
				}
		print OUT "\n";
	}
}

close IN;
close OUT;

open L, ">fusion_len.txt";
foreach (@fusion){
	print L ">$_\n";
	foreach $xx(@{$_}){
		print L "$xx\t";
		}
	print L "\n";
	}
close L;
