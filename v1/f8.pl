use Getopt::Std;
use vars qw( $opt_g $opt_r $opt_t $opt_l );

# Usage
my $usage = "

f8.pl - extract exon sequences that are next to the fusion sites to run blastn.
      by  Qinjie Chu (Version: 1.0  Last Modified: 2017-01-06)

Usage: perl f8.pl [options]
 required:
  -g	genome file.
  -r	refFlat file (simplfied).
  -t	temp7.txt that contains the fusion sites and transcript information.
 optional:
  -l	the length of each exon used to run blastn. [ default = 100bp ]

";

getopts('g:r:t:l:');
die $usage unless $opt_g && $opt_r && $opt_t ;

$genome = $opt_g if $opt_g;
$refflat = $opt_r if $opt_r;
$temp = $opt_t if $opt_t;
$len = $opt_l ? $opt_l : 100;

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

open refflat, "$refflat" or die "please input the refFlat file!\n";
while (<refflat>){
	chomp;
	@refflat=split/\t/,$_;
	#start position is 0-based
	@exonstarts=split/,/,$refflat[9];
	#end position is 1-based
	@exonends=split/,/,$refflat[10];
	foreach $i(0..$#exonstarts){
		$exonstart=$exonstarts[$i]+1;
		$num=$i+1;
		$halocistart{"$refflat[0]\texon$num"}=$exonstart;
		$halociend{"$refflat[0]\texon$num"}=$exonends[$i];
		$hachrom{$refflat[0]}=$refflat[2];
		$haori{$refflat[0]}=$refflat[3];
		}
	}
close refflat;

open temp7, "$temp" or die "please input temp7.txt";
##format of temp7.txt
#LOC_Os01g01450|LOC_Os07g48460	2	3	2reads	J00103:38:H3WYKBBXX:6:1217:24556:14045 1:N:0:ACTTGA,J00103:38:H3WYKBBXX:6:2206:9719:5112 1:N:0:ACTTGA,
open OUT, ">temp8.txt";
while(<temp7>){
	chomp;
	@aa=split/\t/,$_;
	@bb=split/\|/,$aa[0];
	if(exists $ha{"$hachrom{$bb[0]}"}){
		if($haori{"$bb[0]"} eq '+'){
			#左边部分，如果是正向的基因，取外显子后面部分
			$fusionstart=$halociend{"$bb[0]\texon$aa[1]"}-$len;
			$fusionL=substr($ha{"$hachrom{$bb[0]}"},$fusionstart,$len);#substr: position 0-based
			
			#同个基因的ref部分，当前外显子后面部分和下一个外显子前面部分
			$ref1part1=$fusionL;
			$exonnext=$aa[1]+1;
			$ref1part2=substr($ha{"$hachrom{$bb[0]}"},($halocistart{"$bb[0]\texon$exonnext"}-1),$len);
			#$ref1="$ref1part1"."$ref1part2";
			
			}
		if($haori{"$bb[0]"} eq '-'){
			#左边部分，如果是负向的基因，取外显子前面部分，并反向互补
			$fusionL=substr($ha{"$hachrom{$bb[0]}"},($halocistart{"$bb[0]\texon$aa[1]"}-1),$len);
			$fusionL=reverse($fusionL);
			$fusionL=~tr/ACGTacgt/TGCAtgca/;
			
			#同个基因的ref部分，当前外显子前面部分和上一个外显子后面部分，并反向互补
			$ref1part1=$fusionL;
			$exonlast=$aa[1]-1;
			$refstart=$halociend{"$bb[0]\texon$exonlast"}-$len;
			$ref1part2=substr($ha{"$hachrom{$bb[0]}"},$refstart,$len);
			$ref1part2=reverse($ref1part2);
			$ref1part2=~tr/ACGTacgt/TGCAtgca/;
			}
		}
	if(exists $ha{"$hachrom{$bb[1]}"}){
		if($haori{"$bb[1]"} eq '+'){
			#右边部分，如果是正向的基因，取外显子前面部分
			$fusionR=substr($ha{"$hachrom{$bb[1]}"},($halocistart{"$bb[1]\texon$aa[2]"}-1),$len);
			
			#同个基因的ref部分，上一个外显子的后面部分和当前外显子的前面部分
			$exonlast=$aa[2]-1;
			$refstart=$halociend{"$bb[1]\texon$exonlast"}-$len;
			$ref2part1=substr($ha{"$hachrom{$bb[1]}"},$refstart,$len);
			$ref2part2=$fusionR;
			
			}
		if($haori{"$bb[1]"} eq '-'){
			#右边部分，如果是负向的基因，取外显子后面部分，并反向互补
			$fusionstart=$halociend{"$bb[1]\texon$aa[2]"}-$len;
			$fusionR=substr($ha{"$hachrom{$bb[1]}"},$fusionstart,$len);
			$fusionR=reverse($fusionR);
			$fusionR=~tr/ACGTacgt/TGCAtgca/;
			
			#同个基因的ref部分，下一个外显子的前面部分和当前外显子的后面部分，并反向互补
			$exonnext=$aa[2]+1;
			$ref2part1=substr($ha{"$hachrom{$bb[1]}"},($halocistart{"$bb[1]\texon$exonnext"}-1),$len);
			$ref2part1=reverse($ref2part1);
			$ref2part1=~tr/ACGTacgt/TGCAtgca/;
			$ref2part2=$fusionR;
			}
		}
	
		open fusion1, ">fusion1.fa";
		print fusion1 ">fusion1\n$fusionR\n";
	
		open ref1, ">ref1.fa";
		print ref1 ">ref1\n$ref1part2\n";
		
	#blastn
	system("makeblastdb -in ref1.fa -dbtype nucl -out seq_identity");
	system("blastn -db seq_identity -query fusion1.fa -out blastn_out.txt -evalue 1e-5 -outfmt 6 -num_threads 10");

	open BLASTOUT, "blastn_out.txt";
	@lines=<BLASTOUT>;
	@cc=split/\t/,$lines[0];
	print OUT "$aa[0]\t$aa[1]\t$aa[2]\t$cc[2]\t$cc[3]\t";
	@cc=();
	close BLASTOUT;
	
		open fusion2, ">fusion2.fa";
		print fusion2 ">fusion2\n$fusionL\n";
		
		open ref2, ">ref2.fa";
		print ref2 ">ref2\n$ref2part1\n";
	
	#blastn
	system("makeblastdb -in ref2.fa -dbtype nucl -out seq_identity");
	system("blastn -db seq_identity -query fusion2.fa -out blastn_out.txt -evalue 1e-5 -outfmt 6 -num_threads 10");

	open BLASTOUT, "blastn_out.txt";
	@lines=<BLASTOUT>;
	@cc=split/\t/,$lines[0];
	print OUT "$cc[2]\t$cc[3]\t$aa[3]\t$aa[4]\n";
	@cc=();
	close BLASTOUT;
	
	$fusionL="";
	$fusionR="";
	$ref1part1="";
	$ref1part2="";
	$ref2part1="";
	$ref2part2="";
	
	}
close temp7;
close OUT;

	system("rm -f ref1.fa");
	system("rm -f ref2.fa");
	system("rm -f fusion1.fa");
	system("rm -f fusion2.fa");
	system("rm -f blastn_out.txt");
	system("rm -f seq_identity*");