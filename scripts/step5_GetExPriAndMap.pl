#!/usr/bin/perl

# Map the degradome to extended pri-miRNA;

my $usage="perl step5_GetExPriAndMap.pl <proc> <transcriptome> <deg_dir> <deg_list>\n

parameters:
proc          :Number of processors to launch;
transcriptome  :File, transcriptome involved in current task;
deg_dir       :Directory includes the normalized degradome files;
deg_list       :Names of degradome files split by ¡®/¡¯;
\n";

(my $proc,my $tran_file,my $deg_dir,my $deg_list)=@ARGV;
if ($#ARGV!=3) {
	die $usage;
}

my $exadd_dir="miRDP_sort_ex";
my $bwt_dir="bwt_dir";

my $pred_dir;
my $prefix;
if ($tran_file=~/^([^\.]+)\.?/) {
	$prefix=$1;
	$pred_dir=$prefix."_dir";
}


if (!-d "$pred_dir/$exadd_dir" || ! -d $deg_dir) {
	die "Please make sure the directory is there.\n";
}
if (! -d "$pred_dir/$bwt_dir") {
	system(qq(mkdir -p $pred_dir/$bwt_dir));
}


opendir(D,"$pred_dir/$exadd_dir");
my @all_infiles=readdir(D);
closedir(D);

my @infiles;
foreach my $tmpfile (@all_infiles) {
	if ($tmpfile=~/^\./) {
		next;
	}
	else{
		push(@infiles,$tmpfile);
	}
}


my %tag;
foreach my $file (@infiles) {
	open(IN,"$pred_dir/$exadd_dir/$file") || die "step5: $pred_dir/$exadd_dir/$file not found.\n";
	<IN>;
	while (<IN>) {
		chomp;
		my @inf=split;
		$tag{$inf[-1]}=1;
	}
	close(IN);
}

my $expriseq=$prefix."_Extended_pri.fasta";
my $seq_num=1;
my %seqid;
open(O,">$pred_dir/$bwt_dir/$expriseq");
foreach  (keys %tag) {
	my $pri_id="pri-seq_".$seq_num;
	print O ">",$pri_id,"\n",$_,"\n";
	$seqid{$_}=$pri_id;
	$seq_num++;
}
close(O);

###make index
system(qq(bowtie-build $pred_dir/$bwt_dir/$expriseq $pred_dir/$bwt_dir/$expriseq));

my @degfiles=split(/\//,$deg_list);


my %bwtfile;
foreach my $degin (@degfiles) {
	if ($degin=~/^([^\.]+)\.?/) {
		my $mapout=$prefix."_Extended_pri_".$1.".bwt";
		$bwtfile{$mapout}=1;
		system(qq(bowtie -f -v 0 -p $proc -a --norc $pred_dir/$bwt_dir/$expriseq $deg_dir/$degin >$pred_dir/$bwt_dir/$mapout));
	}
}


###BEGIN PARSE
my %lable;
my %match;
foreach  (keys %bwtfile) {
	open(IN,"$pred_dir/$bwt_dir/$_");
	while (<IN>) {
		chomp;
		my @mapresult=split;
		my $map_start=$mapresult[3]+1;
		my $inf=$mapresult[2]."_".$map_start;
		if (!exists $lable{$inf}) {
			$match{$mapresult[2]}.="_".$map_start;
			$lable{$inf}=1;
		}
		else{
			next;
		}

	}
	close(IN);
}

###END PARSE

#------------------add deg Signal---------------
foreach my $file (@infiles){
	my $outfile=$file.".degsignal";
	open(O,">$pred_dir/$bwt_dir/$outfile");
	open(IN,"$pred_dir/$exadd_dir/$file");
	while (<IN>) {
		chomp;
		if (/^miR_mat/) {
			print O $_,"\t","degradation_sigal\n";
		}
		else{
			print O $_,"\t";
			my @subinf=split;
			if (exists $match{$seqid{$subinf[-1]}}) {
				print O $match{$seqid{$subinf[-1]}},"\n";
			}
			else{
				print O "NULL\n";
			}
			
		}
	}
	close(IN);
	close(O);
}
#------------------add deg END------------------


print "step5 finished.\n";