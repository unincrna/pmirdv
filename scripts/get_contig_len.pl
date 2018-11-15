#!/usr/bin/perl

#usage: perl $0 <cdna_dir> <transcriptome>;

my ($trans_dir,$transcriptome)=@ARGV;

my $lenstat;
my %len;
my $id;

if ($transcriptome=~/^([^\.]+)\.?/) {
	$lenstat=$1.".len";
}


open(IN,"$trans_dir/$transcriptome") || die "file $transcriptome cannot be found\n";
while (<IN>) {
	chomp;
	if (/^>(\S+)/) {
		$id=$1;
	}
	else{
		my $sublen=length($_);
		$len{$id}+=$sublen;
	}
}
close(IN);

open(O,">$trans_dir/$lenstat");
foreach  (keys %len) {
	print O $_,"\t",$len{$_},"\n";
}
close(O);