#!/usr/bin/perl

# This script is used to remove the temporary results
# under directories "miRDP","miRDP_parse","miRDP_sort",
# "miRDP_sort_ex","bwt_dir","tmp","pri_exp";

my $usage="perl clean.pl <transcriptome>\n";
if ($#ARGV!=0) {
	die $usage;
}

my $trans_file=$ARGV[0];

my $pred_dir;

if ($trans_file=~/^([^\.]+)\.?/) {
	$pred_dir=$1."_dir";
}
my @dirs=("miRDP","miRDP_parse","miRDP_sort","miRDP_sort_ex","bwt_dir","tmp","pri_exp");
foreach my $dir (@dirs) {
	if (-d "$pred_dir/$dir") {
		system(qq(rm -rf $pred_dir/$dir));
	}
}
