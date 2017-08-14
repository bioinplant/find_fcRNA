find_fcRNA.pl - find fusion-circRNAs (fcRNAs) from RNA-seq data.
      by  Qinjie Chu (Version: 1.0  Last Modified: 2017-04-20 E-mail:qinjiechu@zju.edu.cn)

Usage: perl find_fcRNA.pl [options]
 required:
  -p	path of this package (the end shouldn't be '/').
  -f	text file that contains the names and paths of RNA-Seq data files (fastq format). A name with path of a line in text file. For this text file, it is better to use 'dos2unix' before runing.
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


Installation Prerequisite
1. bowtie 
2. tophat TopHat v2.1.0 
3. perl 
4. samtools

Installation
Unzip the package.

Some tips should be paid attention to:
1. Our 'find_fcRNA' toolkit is suitable for both single-end and paired-end reads.
2. The refFlat file used here have been rearranged. 
   You can find the script in the directory 'refFlat' (refFlat_to_V2.pl), and use it like this: 'perl refFlat_to_V2.pl [old_refFlat.txt]'. Then you can get the V2 version of the refFlat file. It is better to use the V2 version of the refFlat file running 'find_fcRNA'.
   The reason to rearrange the refFlat file is that, only keep one transcript (which have the most exons or longest) for each gene, which could be more convenient and faster for running the scripts.



================================= ATTENTION =========================================
== 1. If you don't want to know the detailed information of running,               ==
==    you don't have to read the following information.                            ==
== 2. You can find the format of each file below.                                  ==
== 3. You can find the final result (fcRNA candidates) in temp8.txt and            ==
==    the parent gene pairs that might create fcRNA candidates in temp7_name.txt . ==
== 4. The sequence similarity of final result in temp8.txt provide users an        ==
==    understanding of these fcRNA candidates. In general, the fcRNA candidates    ==
==    with lower sequence similarity are of higher confidence.                     ==
================================= ATTENTION =========================================


###detailed running steps
1. generate fusion sites by tophat-fusion (v2.1.0) (paired-end reads treated as single-end reads)
   #command£º
      tophat -o [tophat_fusion_output_file] -p [threads] --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search --fusion-min-dist 100000 --fusion-anchor-length 10 --fusion-multireads 10 [genome_index_bowtie1] [one_seq_file]
   #output files:
      [fcRNA_out]/[tophat_fusion_output_file]
                         |----accepted_hits.bam
                         |----deletions.bed
                         |----fusions.out
                         |----insertions.bed
                         |----junctions.bed
                         |----prep_reads.info
                         |----unmapped.bam
                         |----logs(directory)

2. generate unmapped reads by bowtie (v0.12.9)
   #command (if the sequence file is fastq format)£º
      bowtie -t -q -v 2 -p [threads] --mm --un [unmapped_reads_file] [genome_index_bowtie1] [one_fastq_file] > [bowtie_output_file]
   #command (if the sequence file is fasta format)£º
      bowtie -t -f -v 2 -p [threads] --mm --un [unmapped_reads_file] [genome_index_bowtie1] [one_fasta_file] > [bowtie_output_file]
   #output file:
      [fcRNA_out]/[unmapped_reads_file]

3. transfer [tophat_fusion_output_file]/accepted_hits.bam to [tophat_fusion_output_file]/accepted_hits.sam
   #command:
      samtools view -h -o [tophat_fusion_output_file]/accepted_hits.sam [tophat_fusion_output_file]/accepted_hits.bam
   #output file:
      [fcRNA_out]/[tophat_fusion_output_file]/accepted_hits.sam (which will be deleted after f0_simplified_sam.pl)

4. extract useful (fusion) information from accepted_hits.sam file and change the sam file to an simplified one by f0_simplified_sam.pl
   #command:
      perl /path/to/find_fcRNA_package/f0_simplified_sam.pl [tophat_fusion_output_file]/accepted_hits.sam
   #output file:
      [fcRNA_out]/simplified_accepted_hits.sam

5. calculate the exact fusion site of each reads from [fcRNA_out]/simplified_accepted_hits.sam by f1_fusion_site_from_simplified.pl
   #command:
      perl /path/to/find_fcRNA_package/f1_fusion_site_from_simplified.pl -s [fcRNA_out]/simplified_accepted_hits.sam -l [reads_length]
   #usage of f1_fusion_site_from_simplified.pl: 
      perl f1_fusion_site_from_simplified.pl [options]
         required:
            -s	simplified sam file.
            -l	read length.
   #output file:
      [fcRNA_out]/temp1.txt
   #format of temp1.txt
      column 1: chrA-chrB
      column 2: fusion site on chrA
      column 3: fusion site on chrB

6. judge if the fusion sites are inside the transcripts by f2_genic.pl
   #command (fast mode):
      perl /path/to/find_fcRNA_package/f2_genic.pl -r [refFlat_file] -t temp1.txt
   #command (whole mode with '-w'):
      perl /path/to/find_fcRNA_package/f2_genic.pl -r [refFlat_file] -t temp1.txt -w
   #usage of f2_genic.pl
      perl f2_genic.pl [options]
         required:
            -r	refFlat file.
            -t	temp1.txt that contains the fusion sites.
         optional:
            -w	with this option on, while finding fusion sites are genic, it will keep finding others. Default is that, it won't find other gene pairs. 
   #output file:
      [fcRNA_out]/temp2.txt
   #format of temp2.txt
      column 1: chrA-chrB
      column 2: fusion site on chrA
      column 3: the name of the transcript that contains the fusion site on chrA
      column 4: one of the two flanking exons that contains the fusion site on chrA
      column 5: the other of the two flanking exons that contains the fusion site on chrA
      column 6: fusion site on chrB
      column 7: the name of the transcript that contains the fusion site on chrB
      column 8: one of the two flanking exons that contains the fusion site on chrB
      column 9: the other of the two flanking exons that contains the fusion site on chrB

7. extract exon sequences that are adjacent to the fusion sites to run blastn by f3_blastn.pl, and keep the fusion sites with sequence similarity percentage and length less than a centain value.
   #command:
      perl /path/to/find_fcRNA_package/f3_blastn.pl -g [genome_fasta] -r [refFlat_file] -t temp2.txt -l [basepairs] -i [blastn_percentage] -b [blastn_length]
   #usage of f3_blastn.pl
      perl f3_blastn.pl [options]
         required:
            -g	genome file.
            -r	refFlat file (simplfied).
            -t	temp2.txt that contains the fusion sites and transcript information.
         optional:
            -l	the length of each exon used to run blastn. [ default = 100 ]
            -i	identify percentage of blastn result. [ default = 70 ]
            -b	identify length of blastn result. [ default = 100 ]
   #output file:
      [fcRNA_out]/temp3.txt
   #format of temp3.txt
      column 1: the name of the transcript that contains the fusion site on chrA
      column 2: the name of the transcript that contains the fusion site on chrB

8. extract exon sequences of transcript pairs by f4_transcript_exon_seq_rc.pl
   #command:
      perl /path/to/find_fcRNA_package/f4_transcript_exon_seq_rc.pl -g [genome_fasta] -r [refFlat_file] -t temp3.txt
   #usage of f4_transcript_exon_seq_rc.pl
      perl f4_transcript_exon_seq_rc.pl [options]
         required:
            -g	genome file.
            -r	refFlat file (simplfied).
            -t	temp3.txt that contains the fusion transcript pairs.
   #output file:
      [fcRNA_out]/temp4.txt
   #format of temp4.txt (fasta format, the arrangement of exon sequences could be found at http://ibi.zju.edu.cn/bioinplant/)

9. generate exon loci of each sequences in temp4.txt by f5_exon_loci.pl
   #command:
      perl /path/to/find_fcRNA_package/f5_exon_loci.pl -r [refFlat_file]
   #usage of f5_exon_loci.pl
      perl f5_exon_loci.pl [options]
         required:
            -r	refFlat file (simplfied)
   #output file:
      [fcRNA_out]/temp5.txt
      [fcRNA_out]/fusion_len.txt
   #format of temp5.txt
      ##------------------------##------------------------------------------------------------------------------------##
      ##        EXAMPLE         ##             DESCRIPTION                                                            ##
      ##------------------------##------------------------------------------------------------------------------------##
      ##>AAAS|CEP128            ##>transcriptA|transcriptB                                                            ##
      ##258	AAAS-1|CEP128-1     ##junction site of the two exons following	transcript-exonNo.|transcript-exonNo.     ##
      ##------------------------##------------------------------------------------------------------------------------##
   #format of fusion_len.txt
      ##----------------------------------------------------------##
      ##             DESCRIPTION                                  ##
      ##----------------------------------------------------------##
      ##>transcriptA|transcriptB                                  ##
      ##junction site of the every two exons separated by tab     ##
      ##----------------------------------------------------------##

10. build bowtie index of temp4.txt
   #command:
      bowtie-build temp4.txt [bt1_outfile_prefix]
   #output file:
      [fcRNA_out]/[bt1_outfile_prefix].1.ebwt
      [fcRNA_out]/[bt1_outfile_prefix].2.ebwt
      [fcRNA_out]/[bt1_outfile_prefix].3.ebwt
      [fcRNA_out]/[bt1_outfile_prefix].4.ebwt
      [fcRNA_out]/[bt1_outfile_prefix].rev.1.ebwt
      [fcRNA_out]/[bt1_outfile_prefix].rev.2.ebwt

11. map unmapped reads to fusion sequences (temp4.txt) by bowtie
   #command (if the sequence file is fastq format):
      bowtie -t -q -v 2 -p [threads] [bt1_outfile_prefix] [unmapped_reads_file] > [unmapped_to_fusion.txt]
   #command (if the sequence file is fasta format):
      bowtie -t -f -v 2 -p [threads] [bt1_outfile_prefix] [unmapped_reads_file] > [unmapped_to_fusion.txt]
   #output file:
      [fcRNA_out]/[unmapped_to_fusion.txt]
   #format of [unmapped_to_fusion.txt]: bowtie1 out format

12. analyse the bowtie1 result by f6_from_bt1.pl
   #command:
      perl /path/to/find_fcRNA_package/f6_from_bt1.pl -b [unmapped_to_fusion.txt] -t temp5.txt -l [support_length]
   #usage of f6_from_bt1.pl
      perl f6_from_bt1.pl [options]
         required:
            -b	bowtie1 result file.
            -t	temp5.txt that contains loci information of sequences in temp4.txt.
         optional:
            -l	the least length of bases that support one side of junctions. [ default = 10]
   #output file:
      [fcRNA_out]/temp6.txt
   #format of temp6.txt
      ##------------------------------------------------------##
      ##                      EXAMPLE                         ##
      ##------------------------------------------------------##
      ##>DNAJC16|TMEM169                                      ##
      ##4306	4405	DNAJC16-2|TMEM169-1	46_DNAJC16|TMEM169    ##
      ##------------------------------------------------------##

      ##------------------------------------------------------------------------------------##
      ##             DESCRIPTION                                                            ##
      ##------------------------------------------------------------------------------------##
      ##>transcriptA|transcriptB                                                            ##
      ##the format of lines that are not start with '>'                                     ##
      ##   column 1: the start site of reads mapped                                         ##
      ##   column 2: the end site of reads mapped                                           ##
      ##   column 3: the junction information of two transcripts,                           ##
      ##             format is, 'transcript-exonNo.|transcript-exonNo.'                     ##
      ##   column 4: name of the read that support this junction                            ##
      ##------------------------------------------------------------------------------------##

13. further analyse the temp6.txt and identify fusion-circRNA candidates by f7_fusion.pl
   #command:
      perl /path/to/find_fcRNA_package/f7_fusion.pl -r [refFlat_file] -t temp6.txt -s [support_reads_number]
   #usage of f7_fusion.pl
      perl f7_fusion.pl [options]
         required:
            -r	refFlat file.
            -t	temp6.txt that contains bowtie1 map information.
        optional:
            -s	the least number of support reads at each junctions. [ default = 2 ]
   #output file:
      [fcRNA_out]/temp7.txt
      [fcRNA_out]/temp7_name.txt
   #format of temp7.txt
      ##-----------------------------------------------------------------------------------------------------------##
      ##                      EXAMPLE                                                                              ##
      ##-----------------------------------------------------------------------------------------------------------##
      ##ACBD3|CT47A1	4	2	4reads	2288_ACBD3|CT47A10,2293_ACBD3|CT47A10,2296_ACBD3|CT47A10,2298_ACBD3|CT47A10,     ##
      ##CT47A1|ACBD3	2	5	4reads	2313_CT47A10|ACBD3,2312_CT47A10|ACBD3,2321_CT47A10|ACBD3,2322_CT47A10|ACBD3,     ##
      ##-----------------------------------------------------------------------------------------------------------##

      ##------------------------------------------------------------------------------------##
      ##             DESCRIPTION                                                            ##
      ##------------------------------------------------------------------------------------##
      ##for each candidate pair, the format of the first line is:                           ##
      ##   column 1: transcriptA|transcriptB                                                ##
      ##   column 2: the fusion exon No. in transcriptA                                     ##
      ##   column 3: the fusion exon No. in transcriptB                                     ##
      ##   column 4: total number of the supported read of this junction                    ##
      ##   column 5: names of the supported read of this junction                           ##
      ##the format of the second line is:                                                   ##
      ##   column 1: transcriptB|transcriptA                                                ##
      ##   column 2: the fusion exon No. in transcriptB                                     ##
      ##   column 3: the fusion exon No. in transcriptA                                     ##
      ##   column 4: total number of the supported read of this junction                    ##
      ##   column 5: names of the supported read of this junction                           ##
      ##ATTENTION: exons of transcripts described in the first line and the second line     ##
      ##compose a fusion-circRNA candidate.                                                 ##
      ##------------------------------------------------------------------------------------##

   #format of temp7_name.txt (names of parent gene pairs that might generate fcRNAs)
      column 1: transcriptA|transcriptB

14. check if the sequences of fusion sites and the flanking sequences of fusion sites are of high similarity (Blastn e-value < 1e-5) by f8.pl
   #command:
      perl /path/to/find_fcRNA_package/f8.pl -g [genome_fasta] -r [refFlat_file] -t temp7.txt -l [basepairs]
   #usage of f8.pl
      perl f8.pl [options]
         required:
           -g	genome file.
           -r	refFlat file (simplfied).
           -t	temp7.txt that contains the fusion sites and transcript information.
         optional:
           -l	the length of each exon used to run blastn. [ default = 100bp ]
   #output file:
      [fcRNA_out]/temp8.txt
   #format of temp8.txt
      ##---------------------------------------------------------------------------------------------------------------##
      ##                      EXAMPLE                                                                                  ##
      ##---------------------------------------------------------------------------------------------------------------##
      ##ACBD3|CT47A1	4	2					4reads	2288_ACBD3|CT47A10,2293_ACBD3|CT47A10,2296_ACBD3|CT47A10,2298_ACBD3|CT47A10, ##
      ##CT47A1|ACBD3	2	5					4reads	2313_CT47A10|ACBD3,2312_CT47A10|ACBD3,2321_CT47A10|ACBD3,2322_CT47A10|ACBD3, ##
      ##---------------------------------------------------------------------------------------------------------------##

      ##------------------------------------------------------------------------------------##
      ##             DESCRIPTION                                                            ##
      ##------------------------------------------------------------------------------------##
      ##for each candidate pair, the format of the first line is:                           ##
      ##   column 1: transcriptA|transcriptB                                                ##
      ##   column 2: the fusion exon No. in transcriptA                                     ##
      ##   column 3: the fusion exon No. in transcriptB                                     ##
      ##   column 4: the sequence identity percenrage of the sequences of fusion sites and  ##
      ##             the flanking sequences of exon in transcriptA                          ##
      ##   column 5: the sequence identity length of the sequences of fusion sites and      ##
      ##             the flanking sequences of exon in transcriptA                          ##
      ##   column 6: the sequence identity percenrage of the sequences of fusion sites and  ##
      ##             the flanking sequences of exon in transcriptB                          ##
      ##   column 7: the sequence identity length of the sequences of fusion sites and      ##
      ##             the flanking sequences of exon in transcriptB                          ##
      ##   column 8: total number of the supported read of this junction                    ##
      ##   column 9: names of the supported read of this junction                           ##
      ##the format of the second line is:                                                   ##
      ##   column 1: transcriptB|transcriptA                                                ##
      ##   column 2: the fusion exon No. in transcriptB                                     ##
      ##   column 3: the fusion exon No. in transcriptA                                     ##
      ##   column 4: the sequence identity percenrage of the sequences of fusion sites and  ##
      ##             the flanking sequences of exon in transcriptA                          ##
      ##   column 5: the sequence identity length of the sequences of fusion sites and      ##
      ##             the flanking sequences of exon in transcriptA                          ##
      ##   column 6: the sequence identity percenrage of the sequences of fusion sites and  ##
      ##             the flanking sequences of exon in transcriptB                          ##
      ##   column 7: the sequence identity length of the sequences of fusion sites and      ##
      ##             the flanking sequences of exon in transcriptB                          ##
      ##   column 8: total number of the supported read of this junction                    ##
      ##   column 9: names of the supported read of this junction                           ##
      ##ATTENTION: 
      ##1. exons of transcripts described in the first line and the second line             ##
      ##compose a fusion-circRNA candidate.                                                 ##
      ##2. if column 4/5/6/7 is empty, it represents that the sequences of fusion sites and ##
      ##the flanking sequences of fusion sites are of low similarity (Blastn e-value > 1e-5)##
      ##------------------------------------------------------------------------------------##

