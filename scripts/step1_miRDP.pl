#!/usr/bin/perl

my $usage="	perl $0 <proc> <cdna_dir> <transcriptome> <sRNA_dir> <sRNA_file_list>\n";
my ($proc,$trans_dir,$transcriptome,$rdp_file_dir,$rdp_files)=@ARGV;
if ($#ARGV!=4) {
	die $usage;
}

my @srnaFiles=split(/\//,$rdp_files);

my $pred_dir;
my $lenstat;
if ($transcriptome=~/^([^\.]+)\.?/) {
	$pred_dir=$1."_dir";
	$lenstat=$1.".len";
	if (!-d $pred_dir) {
		system(qq(mkdir $pred_dir));
		system(qq(mkdir -p $pred_dir/miRDP));
	}
}

#1
system(qq(bowtie-build $trans_dir/$transcriptome $trans_dir/$transcriptome));

foreach my $rdpfile (@srnaFiles){
	my $data;
	if ($rdpfile=~/^([^\.]+)\.?/) {
		$data=$1;
	}
#2
	if (! -d "$pred_dir/miRDP/$data") {
		system(qq(mkdir -p $pred_dir/miRDP/$data));
	}
	
	my $bwtA=$data.".aln";
	system(qq(bowtie -p $proc -a -v 0 -f $trans_dir/$transcriptome $rdp_file_dir/$rdpfile >$pred_dir/miRDP/$data/$bwtA));
#3
	my $bstA=$data.".bst";
	system(qq(convert_bowtie_to_blast.pl $pred_dir/miRDP/$data/$bwtA $rdp_file_dir/$rdpfile $trans_dir/$transcriptome >$pred_dir/miRDP/$data/$bstA));

#4
	my $preseq=$data."_precursors.fa";
	system(qq(excise_candidate.pl $trans_dir/$transcriptome $pred_dir/miRDP/$data/$bstA 250 >$pred_dir/miRDP/$data/$preseq));

#5
	my $structF=$data."_structures";
	system(qq(cat $pred_dir/miRDP/$data/$preseq | RNAfold --noPS >$pred_dir/miRDP/$data/$structF));

#6
	system(qq(bowtie-build $pred_dir/miRDP/$data/$preseq $pred_dir/miRDP/$data/$preseq));

#7
	my $bwtB=$data."_precursors.aln";
	system(qq(bowtie -a -v 0 -f -p $proc $pred_dir/miRDP/$data/$preseq $rdp_file_dir/$rdpfile >$pred_dir/miRDP/$data/$bwtB));

#8
	my $bstB=$data."_precursors.bst";
	system(qq(convert_bowtie_to_blast.pl $pred_dir/miRDP/$data/$bwtB $rdp_file_dir/$rdpfile $pred_dir/miRDP/$data/$preseq >$pred_dir/miRDP/$data/$bstB));

#9
	my $sigfile=$data."_signatures";
	system(qq(sort +3 -25 $pred_dir/miRDP/$data/$bstB >$pred_dir/miRDP/$data/$sigfile));

#10
	my $predict=$data."_predictions";
	system(qq(miRDP.pl $pred_dir/miRDP/$data/$sigfile $pred_dir/miRDP/$data/$structF >$pred_dir/miRDP/$data/$predict));

#11
#	system(qq(perl /home/hzbio/Softwares/miRDP1.3/get_contig_len.pl $trans_dir/$transcriptome));

#12
	my $nrfile=$data."_nr_predictions";
	my $metP=$data."_filter_P_predictions";
	system(qq(rm_redundant_meet_plant.pl $trans_dir/$lenstat $pred_dir/miRDP/$data/$preseq $pred_dir/miRDP/$data/$predict $pred_dir/miRDP/$data/$nrfile $pred_dir/miRDP/$data/$metP));
	
	print "step1: prediction with $transcriptome and $rdpfile finished\n";
}

