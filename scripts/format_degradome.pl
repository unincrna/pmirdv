#!/usr/bin/perl

my $usage="
perl $0 <input dir> <output dir> <cutoff>
<cutoff> : the length (nt) to cut from the 5'-end of degradome-seq reads\n";

(my $degin_dir,my $degout_dir,my $cutoff)=@ARGV;
if ($#ARGV<2) {
	die "$usage";
}
if (! -d $degin_dir) {
	die "$degin_dir cannot be found\n";
}
else{
	if (! -d $degout_dir) {
		system(qq(mkdir $degout_dir));
	}

	opendir(D,"$degin_dir");
	my @degfiles=readdir(D);
	closedir(D);


	foreach my $indegfile (@degfiles) {
		if ($indegfile=~/^([^\.]+)\.?/) {
			my $prefix=$1;
			my $total;
			my %count;
			my $seq_s;
			my $outdegfile=$prefix.".treat.fasta";
			open(INDEG,"$degin_dir/$indegfile") || die "$indegfile cannot be found\n";
			while (<INDEG>) {
				chomp;
				if (/^>/) {
					next;
				}
				elsif(length($_)<$cutoff){
					next;
				}
				else{
					my $subseq=substr($_,0,$cutoff);
					if ($subseq=~/[^AGCTagct]/) {
						next;
					}
					else{
						$total++;
						$count{$subseq}++;
					}
				}
			}
			close(INDEG);

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
}
