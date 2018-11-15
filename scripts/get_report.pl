#!/usr/bin/perl

#generate the report;
my $usage="perl $0 <organism> <transcriptome> <srna_list> <degradome_list>\n";
my ($org,$transcriptome,$srna_list,$deg_list)=@ARGV;
die "$usage" if ($#ARGV !=3);

my $result_dir="result";
my $prefix;
my $pred_dir;
my $sum_list;
my $priseq_file;
if ($transcriptome=~/^([^\.]+)\.?/) {
	$prefix=$1;
	$pred_dir=$prefix."_dir";
	$sum_list=$org."-".$prefix."-prediction.xls";
	$priseq_file=$prefix."_pri-miRNA.seq";
}

my $count_pri=&count($priseq_file);


my $count_rdp_mature=0;
my %mature_unique;
my $count_rdp_mature_conserved=0;
my %mature_conserved_unique;
my %sigpri;
my $sigA=0;
my $sigB=0;
my $sigC=0;
my $sigD=0;

open(RDP,"$pred_dir/$result_dir/$sum_list")|| die "statistics: NO $sum_list was found.\n";
<RDP>;
while (<RDP>) {
	my @inf=split;
	$count_rdp_mature+=2;
	$mature_unique{$inf[1]}=1;
	$mature_unique{$inf[6]}=1;
	if ($inf[4] ne 'NON') {
		$count_rdp_mature_conserved+=1;
		$mature_conserved_unique{$inf[1]}=1;		
	}
	if ($inf[9] ne 'NON') {
		$count_rdp_mature_conserved+=1;
		$mature_conserved_unique{$inf[6]}=1;		
	}

	my @signals=split(/_/,$inf[-1]);
	my $sigcounts=&sigsum(@signals);
	if ($sigcounts>0) {
		$sigpri{$inf[10]}=1;
	}
	if ($signals[0]) {$sigA++;}
	if ($signals[1]) {$sigB++;}
	if ($signals[2]) {$sigC++;}
	if ($signals[3]) {$sigD++;}

}
close(RDP);
my $count_rdp_mature_unique=keys %mature_unique;
my $count_rdp_mature_conserved_unique=keys %mature_conserved_unique;
my $deg_pri=keys %sigpri;


my $homofile=$org."_conserved_miR.xls";
my $count_conserved=0;
open(HOMO,"$pred_dir/$result_dir/$homofile") || die "statistics: NO $homofile was found.\n";
<HOMO>;
while (<HOMO>) {
	$count_conserved++;
}
close(HOMO);


my $report=$org."_".$prefix."_summary.txt";
open(O,">$pred_dir/$result_dir/$report");
print O "
Organism: $org
Transcriptome : $transcriptome
sRNA-seq      : $srna_list
Degradome-seq : $deg_list

********************
miRDeep-P Result: 
********************
Mature miRNAs (unique sequences): $count_rdp_mature ( $count_rdp_mature_unique )
Conserved mature miRNAs (unique sequences): $count_rdp_mature_conserved ( $count_rdp_mature_conserved_unique )
Pri-miRNAs : $count_pri
Pri-miRNAs with degradome signals: $deg_pri
Degradome-supported cropping sites : site1: $sigA ; site2: $sigB ; site3: $sigC ; site4: $sigD


************************************
Homology search result (miRBase21) 
************************************
Mature miRNAs: $count_conserved

";
close(O);



sub count{
	my $filein=shift;
	my $reads_number=0;
	open(IN,"$pred_dir/$result_dir/$filein") || die "statistics: NO $filein was found.\n";
	while (<IN>) {
		if (/^>/) {
			$reads_number++;
		}
	}
	close(IN);
	return $reads_number;
}


sub sigsum {
	my @sigs=@_;
	my $sig_count;
	foreach $sig (@sigs) {
		$sig_count+=$sig;		
	}
	return $sig_count;
}