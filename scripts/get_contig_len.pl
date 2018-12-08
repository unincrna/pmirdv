#!/usr/bin/perl

use Bio::SeqIO;

my $usage="perl get_contig_len.pl <cdna_dir> <transcriptome>
parameters:
cdna_dir       : directory includes the transcriptome files (FASTA);
transcriptome  : file, transcriptome involved in current task;
\n";

my ($trans_dir,$transcriptome)=@ARGV;

if ($#ARGV!=1) {
	die "$usage\n";
}

if (! -d $trans_dir) {
	die "get_contig_len.pl: Directory $trans_dir cannot be found\n";
}

if (! -e "$trans_dir/$transcriptome") {
	die "get_contig_len.pl: $trans_dir/$transcriptome cannot be found\n";
}

my $lenstat;
my %len;
if ($transcriptome=~/^([^\.]+)\.?/) {
	$lenstat=$1.".len";
	#the output file of length statistics;
}

my $in=new Bio::SeqIO(-file=>"$trans_dir/$transcriptome",-format=>'fasta');
while (my $seq=$in->next_seq()) {
	my $id=$seq->id();
	my $nt=$seq->length();
	$len{$id}=$nt;
}

open(O,">$trans_dir/$lenstat");
foreach  (keys %len) {
	print O $_,"\t",$len{$_},"\n";
}
close(O);
