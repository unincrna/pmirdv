#!/usr/bin/perl

my $usage="usage: perl $0 <extension length> <transcriptome> <cdna_dir>\n\n".
		  "<extension length>   : length of pri-seq extension\n".
		  "<transcriptome file> : transcriptome [FASTA]\n".
		  "<cdna_dir>           : directory contains the transcriptome file\n";

(my $exlen,my $transcriptome, my $cdna_dir)=@ARGV;
die "$usage" if ($#ARGV !=2);

my $sorted_dir="miRDP_sort";
my $exadd_dir="miRDP_sort_ex";

my $pred_dir;
my $prefix;
if ($transcriptome=~/^([^\.]+)\.?/) {
	$prefix=$1;
	$pred_dir=$prefix."_dir";
}


if (! -d "$pred_dir/$sorted_dir") {
	die "the $sorted_dir dirctory is not exist\n";
}
if (! -d "$pred_dir/$exadd_dir") {
	system(qq(mkdir -p $pred_dir/$exadd_dir));
}

my %seq;
my $id;
open(S,"$cdna_dir/$transcriptome") || die "step4: $transcriptome cannot be found\n";
while (<S>) {
	chomp;
	if (/^>(\S+)/) {
		$id=$prefix."_".$1;
	}
	else{
		$seq{$id}.=$_;
	}
}
close(S);

opendir(D,"$pred_dir/$sorted_dir");
my @files=readdir(D);
closedir(D);

foreach my $file (@files) {
	if ($file=~/(.+)\.sort/) {
		my $out=$1.".addex";
		open(O,">$pred_dir/$exadd_dir/$out");
		open(I,"$pred_dir/$sorted_dir/$file") || die "step4: cannot find $file under the directory $sorted_dir\n";
		while (<I>) {
			chomp;
			if (/^miR_mat/) {
				print O $_,"\t","strand_vs_transcripts","\t","Extended_pri-seq","\n";
			}
			else{
				print O $_,"\t";
				my @infs=split;
				if ($infs[11]=~/^(.+)_\d+$/) {
					my $trans_id=$prefix."_".$1;
					if ($seq{$trans_id}=~/$infs[12](.*)$/) {
						print O "+","\t";
						my $candidate_ex=$1;
						if (length($candidate_ex)>=$exlen) {
							my $subplus=substr($candidate_ex,0,$exlen);
							print O $infs[12],$subplus,"\n";
						}
						else{
							print O $infs[12],$candidate_ex,"\n";
						}
					}
					elsif(&revcom($seq{$trans_id})=~/$infs[12](.*)$/){
						print O "-","\t";
						my $candidate_ex2=$1;
						if (length($candidate_ex2)>=$exlen) {
							my $subminus=substr($candidate_ex2,0,$exlen);
							print O $infs[12],$subminus,"\n";
						}
						else{
							print O $infs[12],$candidate_ex2,"\n";
						}
					}
					else{
						print O "please check if the right transcriptome file has been selected\n";
					}
				}

			}
		}
		close(I);
		close(O);
	}
}

print "step4 finished.\n";

sub revcom{
	my $seq=shift;
	$seq=~tr/a-z/A-Z/;
	$seq=~tr/ATGCU/TACGA/;
	$seq=reverse $seq;
}

