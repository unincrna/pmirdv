#!/usr/bin/perl

#clean;

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
