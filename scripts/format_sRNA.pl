#!/usr/bin/perl -w

my $usage=" perl $0 <input dir> <output dir> <minlen> <maxlen>\n";

(my $fastadir, my $rdpdir,my $minlen,my $maxlen)=@ARGV;

if ($#ARGV!=3) {
	die $usage;
}

if (! -d $fastadir) {
	die "$fastadir cannot be found \n";
}
else{
	if (! -d $rdpdir) {
		system(qq(mkdir $rdpdir));
	}

	opendir(FD,"$fastadir");
	my @FAfiles=readdir(FD);
	#closedir(FD);
	foreach my $fafile (@FAfiles) {
		if ($fafile=~/^([^\.]+)\.?/) {
			my $total;
			my %count;
			my $seq_s;
			my $prefix=$1;
			my $rdpfile=$prefix.".rdp.fasta";
			open(I,"$fastadir/$fafile") || die "$fafile cannot be found\n";
			while (<I>) {
				if (/^>/) {
					next;
				}
				else{
					chomp;
					if (/[^AGCTagct]/) {
						next;
					}
					elsif(length($_)<$minlen || length($_)>$maxlen){
						next;
					}
					else{
						$total++;
						$count{$_}++;
					}
				}
			}
			close(I);

			open(O,">$rdpdir/$rdpfile");
			foreach (keys %count) {
				$seq_s++;
				my $rpm=sprintf("%.2f",$count{$_}/$total*1000000);
				my $id=$prefix."_".$seq_s."_".$rpm."_x".$count{$_};
				print O ">",$id,"\n",$_,"\n";
			}
			close(O);
		}
	}
}