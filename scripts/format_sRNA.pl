#!/usr/bin/perl

use Bio::SeqIO;

my $usage="perl $0 <input dir> <output dir> <minlen> <maxlen>\n
<input dir>  : data from small RNA sequecing (FASTA);
<output dir> : formated sequence files (FASTA);
<minlen>     : integer, minimum length of candidate miRNAs;
<maxlen>     : integer, maximum length of candidate miRNAs;
\n";

## Note: Make sure all the files deposited in directory <input>
## are small RNA sequences;


(my $fastadir, my $rdpdir,my $minlen,my $maxlen)=@ARGV;

if ($#ARGV!=3) {
	die $usage;
}

if (! -d $fastadir) {
	die "$fastadir cannot be found \n";
}

if (! -d $rdpdir) {
	system(qq(mkdir $rdpdir));
}

#get the list of input files 

opendir(FD,"$fastadir");
my @FAfiles=readdir(FD);
closedir(FD);

# 

foreach my $fafile (@FAfiles) {
	if ($fafile=~/^([^\.]+)\.?/) {
		my $total;
		my %count;
		my $seq_s;
		my $prefix=$1;
		my $rdpfile=$prefix.".rdp.fasta";
		my $in=new Bio::SeqIO(-file=>"$fastadir/$fafile",-format=>'fasta');
		while (my $seq=$in->next_seq()) {
			my $seq_len=$seq->length();
			my $nt=$seq->seq();
			if ($nt=~/[^ATGCUatgcu]/) {
			##ambiguous bases are not allowed;
				next;
			}
			else{
				if ($seq_len>=$minlen && $seq_len<=$maxlen) {
					$total++;
					$count{$nt}++;
				}
			}
		}

		open(O,">$rdpdir/$rdpfile");
		foreach (keys %count) {
			$seq_s++;
			my $rpm=sprintf("%.2f",$count{$_}/$total*1000000);
			my $id=$prefix."_".$seq_s."_".$rpm."_x".$count{$_};
			##"_x".count is required for miRDP;
			print O ">",$id,"\n",$_,"\n";
		}
		close(O);
	}
}