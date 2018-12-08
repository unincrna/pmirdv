#!/usr/bin/perl

#step7. Full information of a job result.

my $usage="perl step7_Suminfor.pl <organism> <transcriptome>\n
parameters:
organism      :Abbreviation of the organism;
transcriptome  :File, transcriptome involved in current task;
\n";

my ($org,$transcriptome)=@ARGV;
if ($#ARGV!=1) {
	die $usage;
}

my $degadd_dir="bwt_dir";
my $outdir="result";

my $pred_dir;
my $prefix;
my $namefile;
if ($transcriptome=~/^([^\.]+)\.?/) {
	$prefix=$1;
	$pred_dir=$prefix."_dir";
	$namefile=$prefix."_miRDP_information.txt";
}
if (! -d "$pred_dir/$outdir") {
	die $usage;
}
my $outfile=$org."-".$prefix."-prediction.xls";
open(O,">$pred_dir/$outdir/$outfile");
print O "matID_5p\tmatSEQ_5p\tBEGIN_5p\tEND_5p\tHomology_5p\tmatID_3p\tmatSEQ_3p\tBEGIN_3p\tEND_3p\tHomology_3p\tPri_ID\tPre_ID\tDeg_signal\n";


my %first;
my %second;
my %pri_id;
my %pre_id;

open(IN,"$pred_dir/$outdir/$namefile") || die "step7: $namefile cannot be found\n";
while (<IN>) {
	chomp;
	my @names=split(/\t/);
	my $longid=$names[1]."_".$names[3]."_".$names[5]."_".$names[7];
	$first{$longid}=$prefix."_".$names[0];
	$second{$longid}=$prefix."_".$names[2];
	$pri_id{$longid}=$prefix."_".$names[4];
	$pre_id{$longid}=$prefix."_".$names[6];
}
close(IN);


opendir(D,"$pred_dir/$degadd_dir") || die "step7: NO $pred_dir/$degadd_dir was found\n";
my @allfiles=readdir(D);
closedir(D);

my @deg_add_files;
foreach my $tmpfile (@allfiles) {
	if ($tmpfile=~/degsignal$/) {
		push(@deg_add_files,$tmpfile);
	}
}

my %annotation;
my %match;
foreach my $infile (@deg_add_files) {
	open(SIGNAL,"$pred_dir/$degadd_dir/$infile") || die "step7: NO $pred_dir/$degadd_dir/$infile was found\n";
	while (<SIGNAL>) {
		chomp;
		if (/^miR_mat/) {
			next;
		}
		else{
			my @subinfs=split(/\t/);
			my $matureU;
			my $matureD;
			my $homoU;
			my $homoD;
			my ($matUposA,$matUposB,$matDposA,$matDposB);
			if ($subinfs[4] eq 'second') {
				$matureU=$subinfs[10];
				$matureD=$subinfs[7];
				$matUposA=$subinfs[8];
				$matUposB=$subinfs[9];
				$matDposA=$subinfs[5];
				$matDposB=$subinfs[6];
				$homoU=$subinfs[1];
				$homoD=$subinfs[0];
			}
			elsif($subinfs[4] eq 'first'){
				$matureU=$subinfs[7];
				$matureD=$subinfs[10];
				$matUposA=$subinfs[5];
				$matUposB=$subinfs[6];
				$matDposA=$subinfs[8];
				$matDposB=$subinfs[9];
				$homoU=$subinfs[0];
				$homoD=$subinfs[1];

			}
			else{die "step7: Maybe problems in previous steps\n";}
			my $longid=$matureU."_".$matureD."_".$subinfs[12]."_".$subinfs[13];
			if ($annotation{$longid}) {
				next;
			}
			else{
				$annotation{$longid}=$first{$longid}."\t".$matureU."\t".$matUposA."\t".$matUposB."\t".$homoU."\t".$second{$longid}."\t".$matureD."\t".$matDposA."\t".$matDposB."\t".$homoD."\t".$pri_id{$longid}."\t".$pre_id{$longid}."\t";
			}

				my @degsites=split(/_/,$subinfs[-1]);
				my %lable;

				foreach my $site (@degsites) {
					$lable{$site}=1;
				}
				my $phaseA=$matUposA;
				my $phaseB=$matUposB+1;
				my $phaseC=$matDposA;
				my $phaseD=$matDposB+1;

				if ($lable{$phaseA}) {
					$match{$longid}{'A'}=1;
				}
				if ($lable{$phaseB}) {
					$match{$longid}{'B'}=1;
				}
				if ($lable{$phaseC}) {
					$match{$longid}{'C'}=1;
				}
				if ($lable{$phaseD}) {
					$match{$longid}{'D'}=1;
				}
		}
	}
	close(SIGNAL);
}

foreach  (sort {$annotation{$a} cmp $annotation{$b}} keys %annotation) {
	print O $annotation{$_};
	if ($match{$_}{'A'}) {
		print O "1";
	}
	else{
		print O "0";
	}
	if ($match{$_}{'B'}) {
		print O "_1";
	}
	else{
		print O "_0";
	}
	if ($match{$_}{'C'}) {
		print O "_1";
	}
	else{
		print O "_0";
	}
	if ($match{$_}{'D'}) {
		print O "_1";
	}
	else{
		print O "_0";
	}
	print O "\n";
}
close(O);


#clean
`rm $pred_dir/$outdir/$namefile`;

print "step7 finished.\n";