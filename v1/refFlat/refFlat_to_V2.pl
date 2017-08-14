open IN, "$ARGV[0]" or die "please input refFlat file!\n";
$name=(split/\./, $ARGV[0])[0];

@lines=<IN>;
##第一次排序
#先按染色体字符由小到大排列，再按转录起始位置由小到大排列，最后按外显子数目由大到小排列
@line=sort{(split/\t/,$a)[2] cmp (split/\t/,$b)[2] || (split/\t/,$a)[4] <=> (split/\t/,$b)[4] || (split/\t/,$b)[8] <=> (split/\t/,$a)[8]}@lines;

##将转录起始位置一样的转录本只留 外显子数目最多 的一个
open OUT, ">$name\_V2.txt";
foreach (@line){
	@aa=split/\t/,$_;
	if(!exists $ha{"$aa[2]\t$aa[4]\t$aa[5]"}){
		$ha{"$aa[2]\t$aa[4]\t$aa[5]"}=1;
		push @out, $_;
		#print OUT "$_";
		}
	}

##第二次排序
#先按基因名称 字符由小到大排列，再按外显子数目由大到小排列
@out2=sort{(split/\t/,$a)[0] cmp (split/\t/,$b)[0] || (split/\t/,$b)[8] <=> (split/\t/,$a)[8]}@out;

##将基因名称一样的只留 外显子数目最多 的一个
foreach (@out2){
	@bb=split/\t/,$_;
	if(!exists $haname{$bb[0]}){
		$haname{$bb[0]}=1;
		push @out3, $_;
		#print OUT "$_";
		}
	}

##第三次排序
#先按染色体字符由小到大排列，再按转录起始位置由小到大排列
@out4=sort{(split/\t/,$a)[2] cmp (split/\t/,$b)[2] || (split/\t/,$a)[4] <=> (split/\t/,$b)[4]}@out3;
foreach (@out4){
	print OUT "$_";
	}