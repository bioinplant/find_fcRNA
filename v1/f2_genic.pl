use Getopt::Std;
use vars qw( $opt_r $opt_t $opt_w);

# Usage
my $usage = "

f2_genic.pl - tell if the fusion sites are inside the transcript;
               if so, extract exon sequences that are next to the fusion sites to run blastn.
      by  Qinjie Chu (Version: 1.0  Last Modified: 2017-01-11)

Usage: perl f2_genic.pl [options]
 required:
  -r	refFlat file.
  -t	temp1.txt that contains the fusion sites.
 optional:
  -w	with this option on, while finding fusion sites are genic, it will keep finding others. Default is that, it won't find other gene pairs. 

";

getopts('r:t:w');
die $usage unless $opt_r && $opt_t ;

$refflat = $opt_r if $opt_r;
$temp = $opt_t if $opt_t;

open refflat, "$refflat" or die "please input the refFlat file!\n";
while (<refflat>){
	chomp;
	@refflat=split/\t/,$_;
	##format
	#DDX11L1	11873	14409	11873,12612,13220,	12227,12721,14409,
	#gene_name	mRNA_start	mRNA_end	exon_starts	exon_ends
	push @{$refflat[2]}, "$refflat[0]\t$refflat[4]\t$refflat[5]\t$refflat[9]\t$refflat[10]";
	}
close refflat;

open temp1, "$temp" or die "please input the fusion-site (temp1.txt) file!\n";
##format of temp1.txt
#chr1-chr1	104636	578379
#chromA-chromB	fusion-site1	fusion-site2
open OUT, ">temp2.txt";
while (<temp1>){
	chomp;
	@site=split/\t/,$_;
	@chrom=split/-/,$site[0];
	foreach (@{$chrom[0]}){
		@simplerefflat1=split/\t/,$_;
		#fusion-site1 is $site[1], mRNA start position is $simplerefflat1[1]+1(0-based start position), mRNA end position is $simplerefflat1[2]
		if($site[1] > $simplerefflat1[1] && ($site[1] <= $simplerefflat1[2])){
			
			#fusion-site2 is $site[2], mRNA start position is $simplerefflat2[1]+1(0-based start position), mRNA end position is $simplerefflat2[2]
			foreach (@{$chrom[1]}){
				@simplerefflat2=split/\t/,$_;
				if($site[2] > $simplerefflat2[1] && ($site[2] <= $simplerefflat2[2])){
					
					if($simplerefflat1[0] ne $simplerefflat2[0]){
					@exonstarts1=split/,/, $simplerefflat1[3];
					@exonends1=split/,/, $simplerefflat1[4];
					foreach $i(0..$#exonstarts1){
						$exonno=$i+1;
						$exonnoo=$exonno+1;
						#如果fusion-site1在外显子上
						if($site[1]>$exonstarts1[$i] && $site[1]<=$exonends1[$i]){
							if($site[1]*2-$exonstarts1[$i] > $exonends1[$i]){print OUT "$site[0]\t$site[1]\t$simplerefflat1[0]\texon$exonno\texon$exonnoo\t";}
							else{print OUT "$site[0]\t$site[1]\t$simplerefflat1[0]\texon$i\texon$exonno\t";}
							if($opt_w == 0){last;}
							#last;
							}
						#如果fusion-site1在内含子上
						if($site[1]>$exonends1[$i] && $site[1]<=$exonstarts1[$i+1]){
							print OUT "$site[0]\t$site[1]\t$simplerefflat1[0]\texon$exonno\texon$exonnoo\t";
							if($opt_w == 0){last;}
							#last;
							}
						}
					
					@exonstarts2=split/,/, $simplerefflat2[3];
					@exonends2=split/,/, $simplerefflat2[4];
					foreach $j(0..$#exonstarts2){
						$exonno=$j+1;
						$exonnoo=$exonno+1;
						#如果fusion-site2在外显子上
						if($site[2]>$exonstarts2[$j] && $site[2]<=$exonends2[$j]){
							if($site[2]*2-$exonstarts2[$j] > $exonends2[$j]){print OUT "$site[2]\t$simplerefflat2[0]\texon$exonno\texon$exonnoo\n";}
							else{print OUT "$site[2]\t$simplerefflat2[0]\texon$j\texon$exonno\n";}
							if($opt_w == 0){last;}
							#last;
							}
						#如果fusion-site2在内含子上
						if($site[2]>$exonends2[$j] && $site[2]<=$exonstarts2[$j+1]){
							print OUT "$site[2]\t$simplerefflat2[0]\texon$exonno\texon$exonnoo\n";
							if($opt_w == 0){last;}
							#last;
							}
						}
					
					}
					#print OUT "$site[0]\t$site[1]\t$simplerefflat1[0]\t$site[2]\t$simplerefflat2[0]\n";
					if($opt_w == 0){last;}
					#last;
					}
				}
		if($opt_w == 0){last;}
		#last;
			}
		}
	}
close temp1;
close OUT;