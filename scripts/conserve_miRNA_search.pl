#!/usr/bin/perl

my $usage="perl $0 <organism><known_miRNAs><transcriptome><sRNA_dir><srna_list>
parameters:
organism       : Abbreviation of the organism;
known_miRNAs  : Plain text file contains known miRNAs;
transcriptome   : File, transcriptome involved in current task;
sRNA_dir      :Directory includes the normalized sRNA files;
sRNA_file_list  :Names of sRNA files split by ¡®/¡¯;
\n";
my ($org,$plantmiR,$trans_file,$srna_dir,$srna_list)=@ARGV;
if ($#ARGV!=4) {
	die $usage;
}

if (!-d $srna_dir) {
	die "conserve_miRNA_search.pl: NO $sRNA_dir was found.\n";
}

my $result_dir="result";

my $prefix;
my $pred_dir;
if ($trans_file=~/^([^\.]+)\./) {
	$prefix=$1;
	$pred_dir=$prefix."_dir";
}


if (!-d "$pred_dir/$result_dir") {
	die "conserve_miRNA_search.pl: Directory $pred_dir/$result_dir is not found\n";
}


my %tag;
my %datasets;
my %conserved_miR;
my %mirbase;
open(DB,"$plantmiR") || die "conserve_miRNA_search.pl: $plantmiR is not found under current directory";
# $plantmiR is organized as
# temporary_ID	miRNA_sequences	miRbase_ID
# e.g.,
# plantmiR1	AGGGGGCAATCTCACCTCAAC	ata-miR9772b-3p&ata-miR9772a-3p

while (<DB>) {
	chomp;
	my @tmp=split;
	$conserved_miR{$tmp[1]}=$tmp[0];
	$mirbase{$tmp[1]}=$tmp[2];
}
close(DB);


my @srnas=split(/\//,$srna_list);
foreach my $srna (@srnas) {
	my $data;
	my $id;
	if ($srna=~/^([^\.]+)\.?/) {
		$data=$1;
	}
	push(@datasets,$data);
	open(IN,"$srna_dir/$srna") || die "conserve_miRNA_search.pl: sRNA files cannot be found under $sRNA_dir when homology searching.\n";
	while (<IN>) {
		chomp;
		if (/^>(\S+)/) {
			$id=$1;
		}
		else{
			if ($conserved_miR{$_}) {
				my @inf=split(/_/,$id);
				${$data}{$_}=$inf[-2];
				$tag{$_}=1;
			}
		}
	}
	close(IN);
}

my $homology_search=$org."_conserved_miR.xls";
open(OUT,">$pred_dir/$result_dir/$homology_search");
print OUT "miRBase_ID\tMature_seq\t";
print OUT join("\t",@datasets);
print OUT "\n";
foreach  (keys %tag) {
	print OUT $mirbase{$_},"\t",$_;
	for(my $i=0;$i<=$#datasets;$i++) {
		if (${$datasets[$i]}{$_}) {
			print OUT "\t",${$datasets[$i]}{$_};
		}
		else{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close(OUT);

print "homology searching done.\n\n\n";
