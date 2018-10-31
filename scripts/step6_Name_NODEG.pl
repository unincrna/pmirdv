#!/usr/bin/perl

my $usage=" perl $0 <orgnism> <transcriptome>\n";

(my $org,my $transcriptome)=@ARGV;

if ($#ARGV!=1) {
	die $usage;
}


my %tag;
my @files;

my $cdna_data;
my $pred_dir;
if ($transcriptome=~/^([^\.]+)\.?/) {
	$cdna_data=$1;
	$pred_dir=$cdna_data."_dir";
}

my $indir="miRDP_sort";
my $outdir="result";

if (! -d "$pred_dir/$indir") {
	die $!;
}

if (! -d "$pred_dir/$outdir") {
	system(qq(mkdir -p $pred_dir/$outdir));
}


opendir(D,"$pred_dir/$indir") || die $!;
my @allfiles=readdir(D);
closedir(D);

foreach my $item (@allfiles) {
	if ($item=~/\.sort$/) {
		push(@files,$item);
	}
}


foreach my $file (@files) {
		open(I,"$pred_dir/$indir/$file") || "step6: cannot found $file.\n";
		<I>;
		while (<I>) {
			my $mature;
			my $pripre;
			my @inf=split;

			if ($inf[4] eq 'second') {
				$mature=$inf[10]."\t".$inf[7];
			}
			elsif($inf[4] eq 'first'){
				$mature=$inf[7]."\t".$inf[10];
			}
			else{
				print "step6: wrong format of input $pred_dir/$indir/$file .\n";
				}

			$pripre=$inf[12]."\t".$inf[13];
			if (!exists ${$mature}{$pripre}) {
				${$mature}{$pripre}=1;
				$tag{$mature}=1;
			}
		}
		close(I);
}

open(O,">$pred_dir/$outdir/uniqlist.txt");
foreach my $der (keys %tag) {
	foreach my $long (keys %{$der}) {
		print O $der,"\t",$long,"\n";

	}
}
close(O);

############### partI END
my $prin=1;
my $pren=1;
my %pri;
my %pre;
my %ctl;
my %cts;
my $ver;
my $seri;

open(I,"$pred_dir/$outdir/uniqlist.txt");
while (<I>) {
	chomp;
	my @inf=split;
	my $mat=$inf[0]."\t".$inf[1];
	if (!exists $pri{$inf[2]}) {
		$pri{$inf[2]}=$prin;
		$prin++;
	}
	if (!exists $pre{$inf[3]}) {
		$pre{$inf[3]}=$pren;
		$pren++;
	}

	if (!exists ${$inf[2]}{$inf[3]}) {
		$cts{$inf[2]}++;
		${$inf[2]}{$inf[3]}=1;
	}
	my $long=$inf[2]."\t".$inf[3];
	$ctl{$long}++;

}
close(I);

foreach  (keys %ctl) {
	if ($ctl{$_}>1) {
		$ver{$_}=1;
	}
}

foreach  (keys %cts) {
	if ($cts{$_}>1) {
		$seri{$_}="a";
	}
}

my $report_priseq=$cdna_data."_pri-miRNA.seq";
my $report_preseq=$cdna_data."_pre-miRNA.seq";

open(OUTPRI,">$pred_dir/$outdir/$report_priseq");
foreach  (sort {$pri{$a} cmp $pri{$b}}keys %pri) {
	my $tmp_priID=$cdna_data."_".$org."-pri-MIR".$pri{$_};
	print OUTPRI ">",$tmp_priID,"\n",$_,"\n";
}
close(OUTPRI);

open(OUTPRE,">$pred_dir/$outdir/$report_preseq");
foreach  (sort {$pre{$a} cmp $pre{$b}}keys %pre) {
	my $tmp_preID=$cdna_data."_".$org."-pre-MIR".$pre{$_};
	print OUTPRE ">",$tmp_preID,"\n",$_,"\n";

}
close(OUTPRE);


my $report=$cdna_data."_miRDP_information.txt";
my %ver;
my %seri;
my %use_seri;

open(L,"$pred_dir/$outdir/uniqlist.txt");
open(O,">$pred_dir/$outdir/$report");
while (<L>) {
	chomp;
	my @tmp=split;
	my $pripre=$tmp[2]."\t".$tmp[3];
	my ($lastpri,$lastpre,$last5p,$last3p);
	if ($cts{$tmp[2]}==1) {
		if($ctl{$pripre}==1){
			$lastpri=$org."-pri-MIR".$pri{$tmp[2]};
			$lastpre=$org."-pre-MIR".$pre{$tmp[3]};
			$last5p=$org."-miR".$pri{$tmp[2]}."-5p";
			$last3p=$org."-miR".$pri{$tmp[2]}."-3p";
		}
		else{
			$lastpri=$org."-pri-MIR".$pri{$tmp[2]};
			$lastpre=$org."-pre-MIR".$pre{$tmp[3]};
			$last5p=$org."-miR".$pri{$tmp[2]}.".".$ver{$pripre}."-5p";
			$last3p=$org."-miR".$pri{$tmp[2]}.".".$ver{$pripre}."-3p";
			$ver{$pripre}++;
		}
	}
	else{
		if ($ctl{$pripre}==1) {#.123;
			$lastpri=$org."-pri-MIR".$pri{$tmp[2]};
			$lastpre=$org."-pre-MIR".$pre{$tmp[3]};
			$last5p=$org."-miR".$pri{$tmp[2]}.$seri{$tmp[2]}."-5p";
			$last3p=$org."-miR".$pri{$tmp[2]}.$seri{$tmp[2]}."-3p";
			$seri{$tmp[2]}++;
		}
		else{
			if (exists $use_seri{$pripre}) {
			$lastpri=$org."-pri-MIR".$pri{$tmp[2]};
			$lastpre=$org."-pre-MIR".$pre{$tmp[3]};
			$last5p=$org."-miR".$pri{$tmp[2]}.$use_seri{$pripre}.".".$ver{$pripre}."-5p";
			$last3p=$org."-miR".$pri{$tmp[2]}.$use_seri{$pripre}.".".$ver{$pripre}."-3p";
			$ver{$pripre}++;
			}
			else{
			$lastpri=$org."-pri-MIR".$pri{$tmp[2]};
			$lastpre=$org."-pre-MIR".$pre{$tmp[3]};
			$last5p=$org."-miR".$pri{$tmp[2]}.$seri{$tmp[2]}.".".$ver{$pripre}."-5p";
			$last3p=$org."-miR".$pri{$tmp[2]}.$seri{$tmp[2]}.".".$ver{$pripre}."-3p";
			$use_seri{$pripre}=$seri{$tmp[2]};
			$seri{$tmp[2]}++;
			$ver{$pripre}++;
			}
		}

	}

	print O $last5p,"\t",$tmp[0],"\t",$last3p,"\t",$tmp[1],"\t",$lastpri,"\t",$tmp[2],"\t",$lastpre,"\t",$tmp[3],"\n";	
	
}
close(L);
close(O);


#clean
`rm $pred_dir/$outdir/uniqlist.txt`;

print "step6 finished.\n".