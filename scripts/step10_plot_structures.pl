#!/usr/bin/perl

#step10. plot structures of mature miRNAs;

my $usage="perl $0 <org><transcriptome>\n";

my ($org,$trans_file)=@ARGV;
if ($#ARGV!=1) {
	die $usage;
}

my $result_dir="result";
my $prefix;
my $pred_dir;
my $sum_list;
my $priseq_file;

if ($trans_file=~/^([^\.]+)\.?/) {
	$prefix=$1;
	$pred_dir=$prefix."_dir";
	$sum_list=$org."-".$prefix."-prediction.xls";
	$priseq_file=$prefix."_pri-miRNA.seq";
}

if (! -d "$pred_dir/$result_dir/structures") {
	system(qq(mkdir -p "$pred_dir/$result_dir/structures"));
}
if (! -d "$pred_dir/tmp") {
	system(qq(mkdir -p "$pred_dir/tmp"));
}

my %seq;
&getseq($priseq_file);

open(LIST,"$pred_dir/$result_dir/$sum_list") || die "step10. NO $pred_dir/$result_dir/$sum_list was found.\n";
<LIST>;
while (<LIST>) {
	chomp;
	my @full_infs=split;
	my $matID=$full_infs[0];
	my $starID=$full_infs[5];
	my $posFP=$full_infs[2]." ".$full_infs[3];
	my $posTP=$full_infs[7]." ".$full_infs[8];
##id
#	my $priseq_id=$prefix."-".$full_infs[10];
	my $priseq_id=$full_infs[10];

#	my $priseq_id=$full_infs[10];
	my $prifile=$priseq_id.".seq";
	my $priseq_fold=$priseq_id.".fold";
	open(O,">$pred_dir/tmp/$prifile");
	print O $seq{$priseq_id};
	close(O);
	system(qq(RNAfold < $pred_dir/tmp/$prifile > $pred_dir/tmp/$priseq_fold));
	my $plotcom="RNAplot --pre \"".$posFP." 5 GREEN omark ".$posTP." 5 RED omark";
	if ($full_infs[-1]=~/\d_\d_\d_\d/) {
		my @degpos=split(/_/,$full_infs[-1]);
		if ($degpos[0]==1) {
			$plotcom.=" ".$full_infs[2]." cmark";
		}
		if ($degpos[2]==1) {
			$plotcom.=" ".$full_infs[7]." cmark";
		}
		if ($degpos[1]==1) {
			my $posB=$full_infs[3]+1;
			$plotcom.=" ".$posB." cmark";
		}
		if ($degpos[3]==1) {
			my $posB=$full_infs[8]+1;
			$plotcom.=" ".$posB." cmark";
		}
	}
	$plotcom.="\"";
	system(qq($plotcom < $pred_dir/tmp/$priseq_fold));
	my $oldps=$priseq_id."_ss.ps";
	my $psout=$full_infs[0]."_".$full_infs[5].".ps";
	system(qq(mv $oldps $psout));
	system(qq(mv $psout $pred_dir/$result_dir/structures/));
}
close(LIST);

print "step10 finished.\n".


sub getseq{
	my $infile=shift;
	open(INSEQ,"$pred_dir/$result_dir/$infile") || die "step10. NO $pred_dir/$result_dir/$infile was found.\n";
	my $id;
	while (<INSEQ>) {
		if (/^>(\S+)/) {
			$id=$1;
			$seq{$id}=$_;
		}
		else{
			$seq{$id}.=$_;
		}
	}
	close(INSEQ);
}
