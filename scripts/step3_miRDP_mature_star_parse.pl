#!/usr/bin/perl

my $usage="perl $0 <known_miRNAs> <transcriptome>";
die $usage if ($#ARGV!=1);

(my $conserve_set,my $trans_file)=@ARGV;

my $mirdp_dir="miRDP_parse";
my $sorted="miRDP_sort";

my %id;
open(C,"$conserve_set") || die "$conserve_set cannot be found\n";
while (<C>) {
	my @pmiR=split;
	$id{$pmiR[1]}=$pmiR[2];
}
close(C);

my $pred_dir;
if ($trans_file=~/^([^\.]+)\.?/) {
	$pred_dir=$1."_dir";
}



if (-d "$pred_dir/$mirdp_dir") {
	system(qq(ls $pred_dir/$mirdp_dir/*_miRDP.txt >$pred_dir/mirdp.list));
}
else{
	die $!;
}

if (! -d "$pred_dir/$sorted") {
	system(qq(mkdir -p $pred_dir/$sorted));
}


open(L,"$pred_dir/mirdp.list");
while (my $file=<L>) {
	chomp;
	if ($file=~/$mirdp_dir\/(.+)_miRDP/) { #from step2, like SRR1103306_GSM1614023_miRDP.txt
		my $pre=$1;
		my $outfile=$pre.".sort";
		my @tmp=split(/_/,$pre);
		my @RDP;
		open(O,">$pred_dir/$sorted/$outfile");
		open(I,"$file") || die "step3: $file not found\n";
		while (<I>) {
			chomp;
			if (/^mature/) {
				print O "miR_mat","\t","miR_star","\t",$_,"\n";
			}
			else{
				@RDP=split;
				my $miR_mat=exists $id{$RDP[5]} ?  $id{$RDP[5]} : "NON";
				my $miR_star=exists $id{$RDP[8]} ? $id{$RDP[8]} : "NON";
				print O $miR_mat,"\t",$miR_star,"\t";
				print O join("\t",@RDP);
				print O "\n";
			}
		}
		close(I);
		close(O);
	}
}
close(L);

system(qq(rm $pred_dir/mirdp.list));

print "step3 finished.\n";