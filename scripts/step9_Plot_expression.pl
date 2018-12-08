#!/usr/bin/perl


use SVG;

my $usage="
perl step9_Plot_expression.pl <organism> <transcriptome> <result_dir> <srna_dir> <srna_list> <proc>\n
parameters:
organism      :Abbreviation of the organism;
transcriptome  :File, transcriptome involved in current task;
sRNA_dir      :Directory includes the normalized sRNA files;
sRNA_file_list  :Names of sRNA files split by ‘/’;
proc         :Number of processors to launch;
\n";

my ($org,$transcriptome,$srna_dir,$srna_list,$proc)=@ARGV;
if ($#ARGV!=4) {
	die "step9 $usage";
}


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

##short reads mapping and parse

#bowtie-build
my $tmp="pri_exp";
system(qq(mkdir -p $pred_dir/$tmp));
system(qq(cp $pred_dir/$result_dir/$priseq_file $pred_dir/$tmp));
system(qq(bowtie-build $pred_dir/$tmp/$priseq_file $pred_dir/$tmp/$priseq_file));


#map
my %bwtfile;
my @srna_files=split(/\//,$srna_list);
my @groups; #head information;
foreach my $in_srna (@srna_files) {
	if (! -e "$srna_dir/$in_srna") {
		die $usage;
	}

	my $data;
	if ($in_srna=~/^([^\.]+)\.?/) {
		$data=$1;
		push(@groups,$data);
	}
	my $mapout=$data."_mapped.bwt";
	$bwtfile{$mapout}=1;
	system(qq(bowtie -f -p $proc -v 0 -a --norc $pred_dir/$tmp/$priseq_file $srna_dir/$in_srna >$pred_dir/$tmp/$mapout));
}

#parse
my %hash;
foreach my $mapfile (keys %bwtfile) {
	my $srna_set;
	if ($mapfile=~/(.+)_mapped\.bwt/) {
		$srna_set=$1;
	}
	open(I,"$pred_dir/$tmp/$mapfile") || die "step9: NO $pred_dir/$tmp/$mapfile was found.\n";
	while (<I>) {
		chomp;
		my @mapinf=split;
		my @seqinf=split(/_/,$mapinf[0]);
		my $rpm=$seqinf[2];
		my $map_start=$mapinf[3]+1;
		my $map_end=$mapinf[3]+length($mapinf[4]);
		for (my $i=$map_start;$i<=$map_end;$i++) {
			$hash{$srna_set}{$mapinf[2]}{$i}+=$rpm;
		}
	}
	close(I);
}


##pri-seq len BEG
my $pri_id;
my %pri_len;
my @priseqs;
open(SEQ,"$pred_dir/$result_dir/$priseq_file");
while (<SEQ>) {
	chomp;
	if (/^>(\S+)/) {
		$pri_id=$1;
		push(@priseqs,$pri_id);
	}
	else{
		$pri_len{$pri_id}+=length($_);
	}
}
close(SEQ);
##pri-seq len END

##output summary table.

system(qq(mkdir -p $pred_dir/$result_dir/priseq_exp_table));

my %sum;
foreach my $priseq (@priseqs) {
	my $summary_exp=$priseq."_expression.xls";
	open(O,">$pred_dir/$result_dir/priseq_exp_table/$summary_exp");
	print O "position","\t";
	print O join("\t",@groups);
	print O "\t","SUM","\n";
	for (my $j=1;$j<=$pri_len{$priseq};$j++) {
		print O $j;
		my $sumrpm=0;
		foreach my $group (@groups) {
			if (exists $hash{$group}{$priseq}{$j}) {
				print O "\t",$hash{$group}{$priseq}{$j};
				$sumrpm+=$hash{$group}{$priseq}{$j};
			}
			else{
				print O "\t","0";
			}
		}
		print O "\t",$sumrpm,"\n";
		$sum{$priseq}{$j}=$sumrpm;
	}
	close(O);
}


##output summary figures.
system(qq(mkdir -p $pred_dir/$result_dir/priseq_exp_fig));

#get mature regions
my %mature;
my %star;
open(ANNO,"$pred_dir/$result_dir/$sum_list") || die "step9: NO miRNA annotation file was found.\n ";
<ANNO>;
while (<ANNO>) {
	chomp;
	my @annoinf=split;
	my $newpri_id=$annoinf[10];
	if (exists $mature{$newpri_id}) {
		$mature{$newpri_id}.="_".$annoinf[2]."_".$annoinf[3];
	}
	else{
		$mature{$newpri_id}=$annoinf[2]."_".$annoinf[3];
	}

	if (exists $star{$newpri_id}) {
		$star{$newpri_id}.="_".$annoinf[7]."_".$annoinf[8];
	}
	else{
		$star{$newpri_id}=$annoinf[7]."_".$annoinf[8];
	}
}
close(ANNO);


my %region;
foreach my $item (@priseqs) {
	my @mature_sites=split(/_/,$mature{$item});
	my @sort_mature=sort {$a<=>$b} @mature_sites;

	my @star_sites=split(/_/,$star{$item});
	my @sort_star=sort {$a<=>$b} @star_sites;

	$region{$item}{'Mature_start'}=$sort_mature[0];
	$region{$item}{'Mature_end'}=$sort_mature[-1];
	$region{$item}{'Star_start'}=$sort_star[0];
	$region{$item}{'Star_end'}=$sort_star[-1];
	

#plot;
	&plotpri($item);

}

print "step9 finished.\n";


################sub###############

sub plotpri{
	my $width=800;
	my $height=300;
	my $svg=new SVG(width=>$width,height=>$height);
	my $plotpri=shift;
	my @matureA=($region{$plotpri}{'Mature_start'},$region{$plotpri}{'Mature_end'});
	my @matureB=($region{$plotpri}{'Star_start'},$region{$plotpri}{'Star_end'});
	my @values;
	
	for (my $i=1;$i<=$pri_len{$plotpri};$i++) {
		push(@values,$sum{$plotpri}{$i});
	}
	my $widpix=$width*0.9/($#values+1);
	my $maxvalue=&max(@values);
	my $fifth=&unity($maxvalue);
	my $hgtpix=$height*0.9/$fifth/5;

##axis
#y-axis
	my $ybar_x=$width*0.05;#第一个bar的x1；
	my $ybar_y1=$height*0.05;
	my $ybar_y2=$height*0.95;
	$svg->line(x1=>$ybar_x,y1=>$ybar_y1,x2=>$ybar_x,y2=>$ybar_y2,stroke=>'black','stroke-width'=>'2');

#	my $fifth=&unity($maxvalue); # y-values
	my $y_scale_unit=abs($ybar_y2-$ybar_y1)/5;
	my $y_scale_x=$ybar_x+$widpix;

	for (my $i=1;$i<=5;$i++) {
		my $y_scale_y=$ybar_y2-$i*$y_scale_unit;
		my $wrd=$fifth*$i;
		$svg->line(x1=>$ybar_x,y1=>$y_scale_y,x2=>$y_scale_x,y2=>$y_scale_y,stroke=>'black','stroke-width'=>'2');
		my $text_x=$ybar_x+1.5*$widpix;
		$svg->text(x=>$text_x,y=>$y_scale_y,'font-size'=>'12')->cdata($wrd);
	}

#x-axis
	my $xbar_x1=$ybar_x;
	my $xbar_x2=$width;
	my $xbar_y=$ybar_y2;
	$svg->line(x1=>$xbar_x1,y1=>$xbar_y,x2=>$xbar_x2,y2=>$xbar_y,stroke=>'black','stroke-width'=>'2');

##axis


#plot bars;
	my $downy=$height*0.05>10 ? 10:$height*0.05;
	my $mature_color="gold";
	my $other_color="blue";
	for (my $p=0;$p<=$#values;$p++) {
		my $fill;
		my $real_p=$p+1;

		if ($real_p>=$matureA[0] && $real_p<=$matureA[1]) {
			$fill=$mature_color;
		}
		elsif($real_p>=$matureB[0] && $real_p<=$matureB[1]){
			$fill=$mature_color;
		}
		else{
			$fill=$other_color;
		}
		my $box_sx=$ybar_x+$p*$widpix;
		my $box_ex=$ybar_x+$real_p*$widpix;
		my $box_y=$xbar_y-$values[$p]*$hgtpix;

		$svg->path(style => "fill:$fill;stroke:black;stroke-width:0.5;",d => "M $box_sx,$xbar_y L $box_sx,$box_y L $box_ex,$box_y L $box_ex,$xbar_y Z");

		if ($real_p%10==0) {
			my $cx=$box_sx;
			my $cy=$xbar_y+$downy+2;
			$svg->text(x=>$cx,y=>$cy,'font-size'=>'10')->cdata($real_p);
			$svg->line(x1=>($box_sx+$box_ex)/2,y1=>$xbar_y,x2=>($box_sx+$box_ex)/2,y2=>$xbar_y-2,stroke=>'black','stroke-width'=>'2');
		}

		if ($real_p%5==0) {
			$svg->line(x1=>($box_sx+$box_ex)/2,y1=>$xbar_y,x2=>($box_sx+$box_ex)/2,y2=>$xbar_y-1,stroke=>'black','stroke-width'=>'2');
		}
	}

	my $outsvg=$plotpri.".expression.svg";
	open(O,">$pred_dir/$result_dir/priseq_exp_fig/$outsvg");
	print O $svg->xmlify;
	close(O);
}



sub max{
	my @items=@_;
	my $max;
	foreach $item (@items) {
		if ($item>$max) {
			$max=$item;
		}
	}
	return $max;
}

sub unity {
	my $max=shift;
	my $fifth;
	my @candidates=("0","0.01","0.05","0.1","0.5","1","5","10","50","100","500","1000","5000","10000","50000","100000","500000","1000000");
	my $i;

	for ($i=0;$i<=$#candidates;$i++) {
		my $j=$i+1;
		if ($max>$candidates[$i] && $max<=$candidates[$j]) {
			$fifth=$candidates[$j]/5;
		}
	}
	return $fifth;
}