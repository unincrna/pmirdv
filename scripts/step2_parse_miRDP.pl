#!/usr/bin/perl

my $usage= "perl $0 <transcriptome> <sRNA_file_list>\n
parameters:
transcriptome    :File, transcriptome involved in current task;
sRNA_file_list   :Names of sRNA files split by ¡®/¡¯;
\n";

(my $transcriptome,my $srna_list)=@ARGV;
if ($#ARGV!=1) {
	die $usage;
}

my $mirdp_parse="miRDP_parse";
my $trans;
my $pred_dir;
if ($transcriptome=~/^([^\.]+)\.?/) {
	$trans=$1;
	$pred_dir=$trans."_dir";
}

if (! -d "$pred_dir/miRDP") {
	die "step2_parse_miRDP.pl: all the result from miRDeep-P are assumed to be under $pred_dir/miRDP\n";
}

if (! -d "$pred_dir/$mirdp_parse") {
	system(qq(mkdir -p $pred_dir/$mirdp_parse));
}


my @srnas;
my @srnafiles=split(/\//,$srna_list);


for(my $i=0;$i<=$#srnafiles;$i++){
	if($srnafiles[$i]=~/^([^\.]+)\.?/){
		my $data=$1;
		push(@srnas,$data);
	}
}

foreach my $srna (@srnas) {
	if ($srna=~/^\./) {
		next;
	}
	else{
		my $pred=$srna."_predictions";
		my $filt=$srna."_filter_P_predictions";
		my $output=$trans."_".$srna."_miRDP.txt";
		open(O,">$pred_dir/$mirdp_parse/$output");

		my (@id_mat,@id_pri,@seq_pre,@seq_pri,@seq_mat,@seq_star,@arm_mat,@strand_mat,@pos_mat_beg,@pos_mat_end,@pos_star_beg,@pos_star_end,@struct_pre,@struct_pri);
		my $uniqid;
		my %inf;
		open(I,"$pred_dir/miRDP/$srna/$pred") || die "step2: $pred_dir/miRDP/$srna/$pred not found\n";
		while (<I>) {
			chomp;
			if (/^mature_query/) {
				@id_mat=split;
			}
			elsif(/^pri_id/){
				@id_pri=split;
			}
			elsif(/^pre_seq/){
				@seq_pre=split;
			}
			elsif(/^pri_seq/){
				@seq_pri=split;
			}
			elsif(/^mature_seq/){
				@seq_mat=split;
			}
			elsif(/^star_seq/){
				@seq_star=split;
			}
			elsif(/^mature_arm/){
				@arm_mat=split;
			}
			elsif(/^mature_strand/){
				@strand_mat=split;
			}
			elsif(/^mature_beg/){
				@pos_mat_beg=split;
			}
			elsif(/^mature_end/){
				@pos_mat_end=split;
			}
			elsif(/^star_beg/){
				@pos_star_beg=split;
			}
			elsif(/^star_end/){
				@pos_star_end=split;
			}
			elsif(/^pre_struct/){
				@struct_pre=split;
			}
			elsif(/^pri_struct/){
				@struct_pri=split;
			}
			elsif(/^star_struct/){
				$uniqid=$id_pri[1]."_".$id_mat[1];
				$ct{$uniqid}++;
				$inf{$uniqid}=$id_mat[1]."\t".$strand_mat[1]."\t".$arm_mat[1]."\t".$pos_mat_beg[1]."\t".$pos_mat_end[1]."\t".$seq_mat[1]."\t".$pos_star_beg[1]."\t".$pos_star_end[1]."\t".$seq_star[1]."\t".$id_pri[1]."\t".$seq_pri[1]."\t".$seq_pre[1]."\t".$struct_pre[1]."\t".$struct_pri[1]."\n";
			}
		}
		close(I);

		print O "mature_id\tmature_strand\tmature_arm\tmature_beg\tmature_end\tmature_seq\tstar_beg\tstar_end\tstar_seq\tpri_id\tpri_seq\tpre_seq\tpre_structure\tpri_structure\n";

		open(L,"$pred_dir/miRDP/$srna/$filt") || die "step2: $pred_dir/$srna/$filt not found\n";
		while (<L>) {
			my @filt=split;
			my $qryid=$filt[3]."_".$filt[2];
			print O $inf{$qryid};
		}
		close(O);
	}
}

print "step2 finished.\n";