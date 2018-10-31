#!/usr/bin/perl

#step7. Full information of a job result.

my $usage=" perl $0 <organism> <transcriptome>\n";
my ($org,$transcriptome)=@ARGV;
if ($#ARGV!=1) {
	die $usage;
}

my $indir="miRDP_sort";
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
print O "matID_5p\tmatSEQ_5p\tBEGIN_5p\tEND_5p\tHomology_5p\tmatID_3p\tmatSEQ_3p\tBEGIN_3p\tEND_3p\tHomology_3p\tPri_ID\tPre_ID\tDegsignal\n";


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


opendir(D,"$pred_dir/$indir");
my @infiles=readdir(D);
closedir(D);

my %tag; 
foreach my $infile (@infiles) {
		if ($infile=~/^\./) {
			next;
		}
		else{
			my $matposA;my $matposB;
			my $starposA;my $starposB;
			my $matseq;	my $starseq;
			my $priseq; my $preseq;
			my $homomat; my $homostar;
			open(I,"$pred_dir/$indir/$infile") || die "step_7.2: $pred_dir/$indir/$infile not found.\n";
			<I>;
			while (<I>) {
				my @infs=split;
				if ($infs[4] eq 'first') {
					$matposA=$infs[5];
					$matposB=$infs[6];
					$matseq=$infs[7];
					$starposA=$infs[8];
					$starposB=$infs[9];
					$starseq=$infs[10];	
					$homomat=$infs[0];
					$homostar=$infs[1];
				}
				elsif($infs[4] eq 'second'){
					$matposA=$infs[8];
					$matposB=$infs[9];
					$matseq=$infs[10];
					$starposA=$infs[5];
					$starposB=$infs[6];
					$starseq=$infs[7];
					$homomat=$infs[1];
					$homostar=$infs[0];
				}
				$priseq=$infs[12];
				$preseq=$infs[13];

				my $full_id=$matseq."_".$starseq."_".$priseq."_".$preseq;
				if (exists $tag{$full_id}) {
					next;
				}
				else{
					$tag{$full_id}=$first{$full_id}."\t".$matseq."\t".$matposA."\t".$matposB."\t".$homomat."\t".$second{$full_id}."\t".$starseq."\t".$starposA."\t".$starposB."\t".$homostar."\t".$pri_id{$full_id}."\t".$pre_id{$full_id}."\t"."0_0_0_0\n";
				}

			}
			close(I);
		}

}

foreach  (sort {$tag{$a} cmp $tag{$b}} keys %tag) {
	print O $tag{$_};
}

close(O);

#clean
`rm $pred_dir/$outdir/$namefile`;

print "step7 finished.\n".