#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
#GetOptions
# ------------------------------------------------------------------
my ($infile,$ref,$outdir);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$infile,
				"r:s"=>\$ref,
				"o:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($infile );
#######################################################################################
# ------------------------------------------------------------------
# Main Body
# ------------------------------------------------------------------
=pod
mkdir ($outdir) unless (-d $outdir);
mkdir ("$outdir/fasta") unless (-d "$outdir/fasta");
=cut
###############
=pod
my $cmd;
$cmd="formatdb -i $infile -p F -o F";
system($cmd);
$cmd="blastn -query $ref -db $infile -evalue 1e-5 -outfmt 5 -max_hsps 1 -out $outdir/blastn.result.xml -word_size 7";
system($cmd);
$cmd="perl $Bin/blast_parser.pl -tophit 1 -topmatch 1 -m 7 $outdir/blastn.result.xml > $outdir/blastn.tophit.result.xls";
system($cmd);
=cut

###############
my (%cds_subject,%cds_query,%cds_len);
open (IN,"$infile") or die $!;#输入文件
$/="\n"; #定义为换行符
my $id;
while (<IN>) {
	chomp;    #chomp用来去除行尾的$/
	next if (/^$/); #如果没有换行符就跳出循环
	if (/>/) {  #开头有>的话
        #print "$_\n";
		$id=$_;#$_表示默认输出
        #print "$id\n";

	}else{
		$cds_subject{$id}.=$_; #id对应   序列
		#print "$_\n";
		$cds_len{$id} = length $cds_subject{$id};#id--lenth
		#print "$cds_len{$id}\n"
	}
}

open (IN,"$ref") or die $!;
$id=();
while (<IN>) {
	chomp;
	next if (/^$/);
	if (/>/) {
		$id=$_;
	}else{
		$cds_query{$id}.=$_;
		
	}
}

###############
my %cds_homo;
open (IN,"$outdir/blastn.tophit.result.xls") or die $!;
<IN>;#跳过第一行
while (<IN>) {
	chomp;
	my ($qid,$q_length,$q_aln_start,$q_aln_end,$s_length,$s_aln_start,$s_aln_end,$sid)=(split/\t/,$_)[0,1,2,3,5,6,7,15];

=pod
	if($q_aln_start > $q_aln_end or $s_aln_start > $s_aln_end or ($q_aln_end - $q_aln_start + 1) / $q_length < 0.6 or ($s_aln_end - $s_aln_start + 1)/$s_length < 0.6){
		next;#跳过后咋办??????????????????
	}
=cut
	my ($ncbi,$q_pos)=(split/\s+/,$qid);#被任意数量的空格分开 +表示任意数量
	
	my ($q_start,$q_end)=$q_pos=~/(\d+)/g;#全局匹配数字  ???????
	#print $q_start,$q_end;
	#print "\n";
	my $s_pos=(split/\s+/,$sid)[1];
#	my ($s_start,$s_end)=$s_pos=~/(\d+)/g;
	my ($s_start,$s_end)=$s_pos=~/(\d+)/g;
#	my $dis=abs($q_start-$s_start);
	my $dis=length $cds_query{">".$qid};
	#print ">".$qid  与 上面print $id  相同

	push @{$cds_homo{">".$sid}},[$ncbi,$dis,">".$qid];#数组  suject完整信息的哈希, query的名字 长度 query完整带>带基因名的名字
}

###############
my @output;
foreach my $k (keys %cds_homo) {  #$k为subject完整名字带基因名
	my ($NCBI_id,$pos,$gene,$product)=split/\s+/,$k;
	my ($start,$end)=$pos=~/(\d+)/g;#同上疑惑
	($gene)=$gene=~/gene=(.*?)]/;

	my @fasta=("$k\n".$cds_subject{$k});#把MmG1的名字序列放一串字符串里
	my $len = $cds_len{$k};#G1长度
	print "$k len: $len\n";#物种基因情况
	my %filter;
	my @homo_group=@{$cds_homo{$k}};


	@homo_group=sort { 
		$$a[0] cmp $$b[0] 
		or 
		$$b[1]<=>$$a[1]} @homo_group;#排序规则,对Mmg1的每个subject完整信息排序

	my @filter_homo=grep {
		if (!exists $filter{$_->[0]}) {
			$filter{$_->[0]}=1;
			print "$_->[0] $_->[1]\n";#其他物种的基因长度
		}
	}@homo_group;

	for my $homo (@filter_homo) {
		if (!exists $cds_query{$homo->[2]}) {
			print "ERROR: $homo\n";die;
		}
#		next if(length $cds_query{$homo->[2]} < $len * 0.3 or length $cds_query{$homo->[2]} > $len * 1.7);
		push @fasta,"$homo->[2]\n".$cds_query{$homo->[2]};
	}
	my $fasta_str=join("\n",@fasta);
	push @output,[$start,$gene,$fasta_str];
	
}

###############
my $order;
for my $gene (sort {$$a[0]<=>$$b[0]} @output) {
	$order++;
	my $gene_id=$gene->[1];
	open (OUT,">$outdir/fasta/gene$order.$gene_id.fasta");
	print OUT "$gene->[2]\n";
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------


sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact: zhangxq<zhangxq\@genepioneer.com> <627446718\@qq.com>
Description:

Usage:
  Options:
	-i	<infile>	input fasta file, 
	-r	<infile>	ref fasta file,  
	-o	<outdir>	output dir, 
	-h				Help

USAGE
	print $usage;
	exit;
}
