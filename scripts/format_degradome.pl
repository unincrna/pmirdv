#!/usr/bin/perl

use Bio::SeqIO;

my $usage="perl $0 <input dir> <output dir> <cutoff>\n
<cutoff> : the length (nt) of sub-sequences to cut
           from the 5'-end of degradome-seq reads
\n";

(my $degin_dir,my $degout_dir,my $cutoff)=@ARGV;
if ($#ARGV!=2) {
	die "$usage";
}
if (! -d $degin_dir) {
	die "$degin_dir cannot be found\n";
}
if (! -d $degout_dir) {
	system(qq(mkdir $degout_dir));
}

# get the list of degradome-seq files;
# we assume all the files in directory <input dir>
# are FASTA formated degradome-seq files;

opendir(D,"$degin_dir");
my @degfiles=readdir(D);
closedir(D);

#

foreach my $indegfile (@degfiles) {
	if ($indegfile=~/^([^\.]+)\.?/) {
		my $prefix=$1;
		my $total;
		my %count;
		my $seq_s;
		my $outdegfile=$prefix.".treat.fasta";
		my $in=new Bio::SeqIO(-file=>"$degin_dir/$indegfile",-format=>'fasta');
		while (my $seq=$in->next_seq()) {
			my $seq_len=$seq->length();
			my $subseq;
			if ($seq_len < $cutoff) {
				next;
			}
			else{
				$subseq=$seq->subseq(1,$cutoff);
				if ($subseq=~/[^ATGCUatgcu]/) {
					next;
				}
				else{
					$total++;
					$count{$subseq}++;
				}
			}
		}
		
		open(O,">$degout_dir/$outdegfile");
		foreach  (keys %count) {
			$seq_s++;
			my $rpm=sprintf("%.2f",$count{$_}/$total*1000000);
			my $id=$prefix."_".$seq_s."_".$rpm;
			print O ">",$id,"\n",$_,"\n";
		}
		close(O);
	}
}



