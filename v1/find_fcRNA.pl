use warnings;
use strict;

use Getopt::Std;
use vars qw($opt_a $opt_p $opt_f $opt_r $opt_i $opt_l $opt_g $opt_t $opt_b $opt_x $opt_y $opt_s $opt_n $opt_w);

# Usage
my $usage = "

find_fcRNA.pl - find fusion-circRNAs (fcRNAs) from RNA-seq data.
      by  Qinjie Chu (Version: 1.0  Last Modified: 2017-04-20)

Usage: perl find_fcRNA.pl [options]
 required:
  -p	path of this package (the end shouldn't be '/').
  -f	text file that contains the names and paths of RNA-Seq data files (fastq/fasta format). A name with path of a line in text file. For this text file, it is better to use 'dos2unix' before runing.
  -r	refFlat file.
  -i	bowtie1 index of genome file.
  -l	the length of reads.
  -g	genome file (fasta format).
 optional:
  -t	threads to use [ default = 10 ].
  -b	bases of sequences of that are extracted to run blastn. [ default = 100 ]
  -x	identify percentage of first blastn (to delete homologous fusion sites), exceed which will be delete. [ default = 70 ]
  -y	identify length of first blastn (to delete homologous fusion sites), exceed which will be delete. [ default = 100 ]
  -s	the least length of bases that support one side of junctions. [ default = 10]
  -n	the least number of support reads at each junctions. [ default = 2 ]
  -a	with -a option, it indicates RNA-Seq data are fasta format. [ default: off ]
  -w	with this option on, while finding fusion sites are genic, it will keep finding others. [ default: off (it won't find other gene pairs) ] 

";

my ($path, $file, $refFlat, $genome_index, $len, $genome_fasta);
my ($thread, $basepairs, $blastn_percentage, $blastn_length, $support_length, $support_number, $exon_extract);
my ($opt_aa, $opt_ww);

getopts('ap:f:r:i:l:g:t:b:x:y:s:n:w');
die $usage unless $opt_p && $opt_f && $opt_r && $opt_i && $opt_l && $opt_g;
$path = $opt_p if $opt_p;
$file = $opt_f if $opt_f;
$refFlat = $opt_r if $opt_r;
$genome_index = $opt_i if $opt_i;
$len = $opt_l if $opt_l;
$genome_fasta = $opt_g if $opt_g;

$thread = $opt_t ? $opt_t : 10;
$basepairs = $opt_b ? $opt_b : 100;
$blastn_percentage = $opt_x ? $opt_x : 70;
$blastn_length = $opt_y ? $opt_y : 100;
$support_length = $opt_s ? $opt_s : 10;
$support_number = $opt_n ? $opt_n : 2;
$opt_aa = $opt_a ? $opt_a : 0;
$opt_ww = $opt_w ? $opt_w : 0;

my (@runs, $i, $tophatfusion, $bowtie, $viewbam, $step0, $step1, $step2, $step3, $step4, $step5, $step6, $step7, $step8);
my (@name, $name, $outf);
open IN, "$file" or die "Please input text file that contains the names and paths of RNA-Seq data files (fastq format)! If data files are fasta format, please use the '-a' option.\n";
while(<IN>){
	chomp;
	push @runs, $_;
	}
close IN;


foreach $i(@runs){
	
	# Manage outfile name.
	@name = split(/\//, $i);
	$name = $name[$#name];
	@name = split(/\./, $name);
	$outf = $name[0];
	
	#create out file
	system(qq(mkdir $outf\_fcRNA_out));
	my $wd1 = $ENV{'PWD'};
	chdir "$wd1/$outf\_fcRNA_out/";
	print "The working dir was changed to: $wd1/$outf\_fcRNA_out/\n\n";

	
	##tophat-fusion
	#paired-end reads treated as single-end reads
	$tophatfusion=qq(tophat -o $outf\_tophat_fusion -p $thread --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search --fusion-min-dist 100000 --fusion-anchor-length 10 --fusion-multireads 10 $genome_index $i);
	print "Tophat-fusion is running...\n	The command is: $tophatfusion\n";
	system($tophatfusion);
	print "Tophat-fusion...Done!\n\n";
	
	##bowtie to generate unmapped reads
	if($opt_aa == 0){$bowtie=qq(bowtie -t -q -v 2 -p $thread --mm --un $outf\_unmapped $genome_index $i > $outf\_to_genome.txt );}
	if($opt_aa == 1){$bowtie=qq(bowtie -t -f -v 2 -p $thread --mm --un $outf\_unmapped $genome_index $i > $outf\_to_genome.txt );}
	print "Bowtie to generate reads that are unmapped to the genome, running...\n	The command is: $bowtie\n";
	system($bowtie);
	print "Bowtie to generate unmapped reads...Done!\n\n";
	system(qq(rm -f $outf\_to_genome.txt));
	
	##transfer bam to sam
	$viewbam=qq(samtools view -h -o $outf\_tophat_fusion/accepted_hits.sam $outf\_tophat_fusion/accepted_hits.bam );
	print "Bam file is changing to sam file...\n	The command is: $viewbam\n";
	system($viewbam);
	print "Bam file to sam file... Done!\n\n";
	
	##extract useful information of accepted_hits.sam file and change the sam file to an simplified one
	$step0=qq(perl $path/f0_simplified_sam.pl $outf\_tophat_fusion/accepted_hits.sam);
	print "$step0\n";
	system($step0);
	system(qq(rm -f $outf\_tophat_fusion/accepted_hits.sam));
	
	##start the method
	$step1=qq(perl $path/f1_fusion_site_from_simplified.pl -s simplified_accepted_hits.sam -l $len);
	print "$step1\n";
	system($step1);
	
	if($opt_ww == 0){$step2=qq(perl $path/f2_genic.pl -r $refFlat -t temp1.txt);}
	if($opt_ww == 1){$step2=qq(perl $path/f2_genic.pl -r $refFlat -t temp1.txt -w);}
	print "$step2\n";
	system($step2);
	
	$step3=qq(perl $path/f3_blastn.pl -g $genome_fasta -r $refFlat -t temp2.txt -l $basepairs -i $blastn_percentage -b $blastn_length);
	print "$step3\n";
	system($step3);
	
	$step4=qq(perl $path/f4_transcript_exon_seq_rc.pl -g $genome_fasta -r $refFlat -t temp3.txt );
	print "$step4\n";
	system($step4);
	
	$step5=qq(perl $path/f5_exon_loci.pl -r $refFlat);
	print "$step5\n";
	system($step5);
	
	#build bowtie index
	system(qq(bowtie-build temp4.txt $outf\_bt1));
	
	#unmapped reads map to fusion sequences by bowtie
	if($opt_aa == 0){
		system(qq(bowtie -t -q -v 2 -p $thread $outf\_bt1 $outf\_unmapped > $outf\_unmapped_to_fusion.txt ));
		}
	if($opt_aa == 1){
		system(qq(bowtie -t -f -v 2 -p $thread $outf\_bt1 $outf\_unmapped > $outf\_unmapped_to_fusion.txt ));
		}
	
	$step6=qq(perl $path/f6_from_bt1.pl -b $outf\_unmapped_to_fusion.txt -t temp5.txt -l $support_length );
	print "$step6\n";
	system($step6);
	
	$step7=qq(perl $path/f7_fusion.pl -r $refFlat -t temp6.txt -s $support_number );
	print "$step7\n";
	system($step7);
	
	$step8=qq(perl $path/f8.pl -g $genome_fasta -r $refFlat -t temp7.txt -l $basepairs );
	print "$step8\n";
	system($step8);
	
	system(qq(mv temp6.txt $outf\_temp6.txt));
	system(qq(mv temp7.txt $outf\_temp7.txt));
	system(qq(mv temp7_name.txt $outf\_temp7_name.txt));
	system(qq(mv temp8.txt $outf\_temp8.txt));

	
	#go back to the upper folder
	my $wd2=$ENV{'PWD'};
	$wd2 =~ s/$outf\_fcRNA_out//;
	chdir "$wd2";
	print "The working dir was changed to: $wd2\nFusion-circRNA identification of $name, done!\n\n";


	}




