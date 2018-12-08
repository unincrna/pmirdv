#!/usr/bin/perl


#############################  USAGE  #######################################
my $usage=
"usage: $0 --joblist joblist.txt [options] 

***************This is the main script of PmirDV pipeline.***************** 

PARAMETERS:

--help|-h : Usage
--org     : Organism(abbreviation,like ath);
--joblist : Job information,each task per line [String];
--proc|-p : Number of processors to use [Integer,default 3];
--exlen   : Length of pri-miRNA extension [Integer,default 50];
--cdna    : Directory contains transcriptome files [Default: ./];
--srna    : Directory contains sRNA files [Default: ./];
--deg     : Directory contains degradome files [Default: ./];
--ref     : known miRNAs[Default: plantmiR.list];
";

###########################  USAGE-END  #####################################

#--------------------------getopt--------------------------
use Getopt::Long;
my ($org,$joblist,$transcriptome_dir,$srna_dir,$degradome_dir,$proc,$description,$known_miRNAs);

GetOptions(
	'org:s' =>\$org,
	'joblist=s'=>\$joblist,
	'exlen:i'=>\$exlen,
	'cdna:s'=>\$transcriptome_dir,     
	'srna:s'=>\$srna_dir,              
	'deg:s'=>\$degradome_dir,
	'proc|p:i'=>\$proc,
	'ref=s'=>\$known_miRNAs,
	'help|h'=>\$description,
);

if($description){die $usage;}
if(! $org){die $usage;}
if (! $joblist) {die $usage;}
if (! $transcriptome_dir) {$transcriptome_dir="./";}
if (! $srna_dir) {$srna_dir="./";}
if (! $degradome_dir) {$degradome_dir="./";}
if (! $proc) {$proc=3;}
if (! $exlen) {$exlen=50;}
if (! $known_miRNAs){$known_miRNAs="plantmiR.list";}

#-------------------------getopt-END----------------------


#my $known_miRNAs="plantmiR.list";
if (! -e $known_miRNAs) {
	die "File of the reference miRNAs cannot be found\n";
}

open(JOB,"$joblist") || die "$joblist cannot be found.\n";
while (<JOB>) {
	chomp;
	my @set=split; 
	my $trans_file=$set[0];

	if (! -e "$transcriptome_dir/$trans_file") {
		die "$trans_file cannot be found.\n"; #transcriptome file;
	}	

	&check_file($srna_dir,$set[1]); #sRNA files;
	&check_id($srna_dir,$set[1]);#sRNA files;

	if ($set[2] eq 'NULL') {
		print "No Degradome information provided.\n";
	}
	else{
		&check_file($degradome_dir,$set[2]); #degradome files
	}


#get transcripts length;
	system(qq(perl ./scripts/get_contig_len.pl $transcriptome_dir $trans_file));


#miRDP and parse
	system(qq(perl ./scripts/step1_miRDP.pl $proc $transcriptome_dir $trans_file $srna_dir $set[1]));
	system(qq(perl ./scripts/step2_parse_miRDP.pl $trans_file $set[1]));
	system(qq(perl ./scripts/step3_miRDP_mature_star_parse.pl $known_miRNAs $trans_file));
#miRDP END


#degradation signal
if ($set[2] eq 'NULL') {
	system(qq(perl ./scripts/step6_Name_NODEG.pl $org $trans_file));
	system(qq(perl ./scripts/step7_Suminfor_NODEG.pl $org $trans_file)); 
}
else{
	system(qq(perl ./scripts/step4_add_ex.pl $exlen $trans_file $transcriptome_dir));
	system(qq(perl ./scripts/step5_GetExPriAndMap.pl $proc $trans_file $degradome_dir $set[2]));
	system(qq(perl ./scripts/step6_Name.pl $org $trans_file));
	system(qq(perl ./scripts/step7_Suminfor.pl $org $trans_file)); 
}

	system(qq(perl ./scripts/step8_GetExpressionProfile.pl $org $trans_file $srna_dir $set[1]));
	system(qq(perl ./scripts/step9_Plot_expression.pl $org $trans_file $srna_dir $set[1] $proc));
	system(qq(perl ./scripts/step10_plot_structures.pl $org $trans_file));

#degradation signal  END


#homology search
	system(qq(perl ./scripts/conserve_miRNA_search.pl $org $known_miRNAs $trans_file $srna_dir $set[1]));
#homology search END

#report
	system(qq(perl ./scripts/get_report.pl $org $trans_file $set[1] $set[2]));

#report END


#clean
	system(qq(perl ./scripts/clean.pl $trans_file));

#clean END

}



###sub###

sub check_file {
	my $dir=shift;
	if ($dir=~/(.+)\/$/) {
		$dir=$1;
	}
	my @files=split(/\//,shift);
	foreach my $file (@files) {
		if (-e "$dir/$file") {
			print "$file","\." x 6, "found\n";
		}
		else{
			die "$file cannot be found.\n";
		}
	}
}


sub check_id {
	my $dir=shift;
	my @files=split(/\//,shift);
	foreach my $file (@files) {
		my $flag;
		open(CKID,"$dir/$file") || die "$file cannot be found for ID checking.\n";
		while (<CKID>) {
		if (/^>(\S+)/) {
			my @longid=split(/_/,$1);
			if ($#longid!=3) {
				die "Incorrect format of sequence IDs in $file.\n";
			}
			elsif($longid[3]!~/x\d+/){
				die "Incorrect format of sequence IDs in $file.\n";
			}
			elsif($longid[2]=~/[^\.0-9]/ || $longid[4]=~/[^0-9]/){
				die "Incorrect format of sequence IDs in $file.\n";
			}
			else{$flag=1;}
		}
	  }
	close(CKID);
	if ($flag) {
		print $file,"\." x 6,"checked.\n";
	    }
	}
}

