#!/usr/bin/perl

my $usage="perl $0 <organism><transcriptome><sRNA_dir><sRNA_file_list>\n";

my ($org,$transcriptome,$srna_dir,$srna_list)=@ARGV;

if ($#ARGV!=3) {
	die "step8: $usage";
}


my $result_dir="result";
my $prefix;
my $pred_dir;
my $result;
if ($transcriptome=~/^([^\.]+)\.?/) {
	$prefix=$1;
	$pred_dir=$prefix."_dir";
	$result=$org."-".$prefix."-prediction.xls";
}

if (! -d $pred_dir || ! -d $srna_dir) {
	die "step8: No work directory is already exists. \n";
}

my @srna_files=split(/\//,$srna_list);
my %express;
foreach $srna_data (@srna_files) {
	&GetExp($srna_data);
}

my $exp_profile=$org."-".$prefix."-miRNA_expression.xls";

open(OUT,">$pred_dir/$result_dir/$exp_profile");
print OUT "Mature_ID","\t","Mature_seq","\t","Conservation";
foreach my $insrna (@srna_files) {
	if ($insrna=~/^([^\.]+)\.?/) {
		print OUT "\t",$1;
	}
}
print OUT "\n";

open(LIST,"$pred_dir/$result_dir/$result") || die "step8: $pred_dir/$result_dir/$result cannot be found.\n";
<LIST>;
while (<LIST>) {
	chomp;
	my @record=split;
#5p
	print OUT $record[0],"\t",$record[1],"\t",$record[4],"\t";
	my $dataset;
	my %sum_exp;
	foreach my $item (@srna_files) {
		if ($item=~/^([^\.]+)\.?/) {
			$dataset=$1;
			if (exists $hash{$dataset}{$record[1]}) {
				print OUT $hash{$dataset}{$record[1]},"\t";
				$sum_exp{$record[1]}+=$hash{$dataset}{$record[1]};
			}
			else{
				print OUT "0","\t";
			}
		}
	}
	print OUT "\n";

#3p

	print OUT $record[5],"\t",$record[6],"\t",$record[9],"\t";
	my $dataset;
	foreach my $item (@srna_files) {
		if ($item=~/^([^\.]+)\.?/) {
			$dataset=$1;
			if (exists $hash{$dataset}{$record[6]}) {
				print OUT $hash{$dataset}{$record[6]},"\t";
			}
			else{
				print OUT "0","\t";
			}
		}
	}
	print OUT "\n";
}
close(LIST);

close(OUT);

#

sub GetExp{
	my $infile=shift;
	my $data;
	if ($infile=~/^([^\.]+)\.?/) {
		$data=$1;
	}
	my @srna_inf;
	open(IN,"$srna_dir/$infile") || die "step8: $srna_dir/$infile cannot be found.\n";
	while (<IN>) {
		chomp;
		if (/^>/) {
			@srna_inf=split(/_/);
		}
		else{
			$hash{$data}{$_}=$srna_inf[-2];
		}
	}
	close(IN);
}

print "step8 finished.\n".