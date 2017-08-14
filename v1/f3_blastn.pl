use Getopt::Std;
use vars qw( $opt_g $opt_r $opt_t $opt_l $opt_i $opt_b);

# Usage
my $usage = "

f3_blastn.pl - extract exon sequences that are next to the fusion sites to run blastn.
      by  Qinjie Chu (Version: 1.0  Last Modified: 2017-01-10)

Usage: perl f3_blastn.pl [options]
 required:
  -g	genome file.
  -r	refFlat file (simplfied).
  -t	temp2.txt that contains the fusion sites and transcript information.
 optional:
  -l	the length of each exon used to run blastn. [ default = 100 ]
  -i	identify percentage of blastn result. [ default = 70 ]
  -b	identify length of blastn result. [ default = 100 ]

";

getopts('g:r:t:l:i:b:');
die $usage unless $opt_g && $opt_r && $opt_t ;

$genome = $opt_g if $opt_g;
$refflat = $opt_r if $opt_r;
$temp = $opt_t if $opt_t;
$len = $opt_l ? $opt_l : 100;
$blastnidentify = $opt_i ? $opt_i : 70;
$blastnlen = $opt_b ? $opt_b : 100;

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
		}
	}
close refflat;

open temp2, "$temp" or die "please input temp2.txt";
##format of temp2.txt
#chr1-chr1	139571	LOC729737	exon1	exon2	224139867	GTF2IP20	exon1	exon2
#chromA-chromB	fusion-site1	transcript1	fusion-site1 flanking left exon	fusion-site1 flanking right exon	fusion-site2	transcript2	fusion-site2 flanking left exon	fusion-site2 flanking right exon
while(<temp2>){
	chomp;
	@location=split/\t/,$_;
	push @pairs, "$location[2]\t$location[3]\t$location[4]\t$location[6]\t$location[7]\t$location[8]";
	}
close temp2;

@uniquePairs=grep{++$count{$_}<2}@pairs;
$numUniquePairs=$#uniquePairs+1;
##format of @uniquePairs
#LOC729737	exon1	exon2	GTF2IP20	exon1	exon2

open OUT, ">temp3.txt";
foreach (@uniquePairs){
	@pp=split/\t/,$_;
	
	##extract left fusion site flanking exon sequences
	#$hachrom{$pp[0]} is the chromsomose of left fusion site;
	if(exists $ha{"$hachrom{$pp[0]}"}){
		if($pp[1] eq exon0){$leftL="";}
		if($pp[1] ne exon0){
			$leftstart=$halociend{"$pp[0]\t$pp[1]"}-$len;
			$leftL=substr($ha{"$hachrom{$pp[0]}"},$leftstart,$len);#substr: position 0-based
			}
		$leftR=substr($ha{"$hachrom{$pp[0]}"},$halocistart{"$pp[0]\t$pp[2]"},$len);
		
		open left, ">left.fa";
		print left ">left\n$leftL$leftR\n";
		
		}
	
	##extract right fusion site flanking exon sequences
	# $hachrom{$pp[3]} is the chromsomose of right fusion site
	if(exists $ha{"$hachrom{$pp[3]}"}){
		if($pp[4] eq exon0){$rightL="";}
		if($pp[4] ne exon0){
			$rightstart=$halociend{"$pp[3]\t$pp[4]"}-$len;
			$rightL=substr($ha{"$hachrom{$pp[3]}"},$rightstart,$len);
			}
		$rightR=substr($ha{"$hachrom{$pp[3]}"},$halocistart{"$pp[3]\t$pp[5]"},$len);
		
		open right, ">right.fa";
		print right ">right\n$rightL$rightR\n";
		
		}
	
	##print progressing
	$k++;
	print "This is the $k , total $numUniquePairs pairs.\n";
	
	#blastn
	system("makeblastdb -in left.fa -dbtype nucl -out seq_identity");
	system("blastn -db seq_identity -query right.fa -out blastn_out.txt -evalue 1e-5 -outfmt 6 -num_threads 10");
	
	open BLASTOUT, "blastn_out.txt";
	@lines=<BLASTOUT>;
	@cc=split/\t/,$lines[0];
	#defult: identity less than 70% and length less than 101bp
	if(($cc[2] < $blastnidentify) && ($cc[3] <= $blastnlen)){
		print OUT "$pp[0]\t$pp[3]\n";##fusion transcript pairs
		}
	close BLASTOUT;
	
	}
	

close OUT;

	system("rm -f left.fa");
	system("rm -f right.fa");
	system("rm -f blastn_out.txt");
	system("rm -f seq_identity*");