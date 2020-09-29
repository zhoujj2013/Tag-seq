#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);
use List::Util qw/sum/;

&usage if @ARGV<1;

sub usage {
        my $usage = << "USAGE";

        This script detect DNA editing site by Tag-seq.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 <t_site_plus_plus_f> <t_site_plus_minus_f> <t_site_minus_plus_f> <t_site_minus_minus_f> xx.control.bed grna_f pam blacklist_f prefix
	
	Example:
	t_site_plus_plus_f: t_site_plus_plus.proximal.sorted

USAGE
print "$usage";
exit(1);
};

############## Get input parameters ##############
my $t_site_plus_plus_f=shift;
my $t_site_plus_minus_f=shift;
my $t_site_minus_plus_f=shift;
my $t_site_minus_minus_f=shift;

my $bg_f = shift;

my $grna_f=shift;
my $pam = shift;
my $blacklist_f=shift;

my $ref = shift;
my $chrom_size = shift;

my $prefix=shift;
my $max_mismatch = shift;
my $max_gap = shift;
############## parameters end ####################


##################################################
# merge nearby editting sites
my %tcount;
foreach my $f (($t_site_plus_plus_f,$t_site_plus_minus_f,$t_site_minus_plus_f,$t_site_minus_minus_f)){
open IN,"$f" || die $!;
	while(<IN>){
		chomp;
		my @t = split /\t/;
		$tcount{$t[3]} = $t[-1];
	}
	close IN;
}
#print Dumper(\%tcount);

`cat $t_site_plus_plus_f $t_site_plus_minus_f $t_site_minus_plus_f $t_site_minus_minus_f | bedtools sort -i - > $prefix.all.sites.bed`;
# this is an important cutoff, I will test 5bp and 10bp
`bedtools merge -d 10 -c 4 -o collapse -i $prefix.all.sites.bed > $prefix.all.sites.merged`;

open OUT,">","$prefix.all.sites.merged.confirmed.raw" || die $!;
open IN,"$prefix.all.sites.merged" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $start = int(($t[1]+$t[2])/2);
	my $end = $start + 1;
	my @id = split /,/,$t[3];
	my %c;
	# initiation
	$c{"plus_plus"} = 0;
	$c{'plus_minus'} = 0;
	$c{'minus_plus'} = 0;
	$c{'minus_minus'} = 0;

	foreach my $id (@id){
		my @tt = split /_/,$id;
		my $index = $tt[1]."_".$tt[2];
		$c{$index} = $c{$index} + $tcount{$id};
	}

	# this is the filtering scheme (old)
	#if(($c{"plus_plus"} >= 5 and $c{'plus_minus'} >= 5) or ($c{'minus_plus'} >= 5 and $c{'minus_minus'} >= 5)){
	#	print OUT "$t[0]\t$start\t$end\t$t[3]\t$c{'plus_plus'}\t$c{'plus_minus'}\t$c{'minus_plus'}\t$c{'minus_minus'}\n";
	#}elsif(($c{"plus_plus"} >= 5 or $c{'plus_minus'} >= 5) and ($c{'minus_plus'} >= 5 or $c{'minus_minus'} >= 5)){
	#	print OUT "$t[0]\t$start\t$end\t$t[3]\t$c{'plus_plus'}\t$c{'plus_minus'}\t$c{'minus_plus'}\t$c{'minus_minus'}\n";
	#}
	# this is the filtering scheme (new)
	my @count = ($c{"plus_plus"}, $c{'plus_minus'}, $c{'minus_plus'}, $c{'minus_minus'});
	my $flag = 0;
	foreach my $count (@count){
		if($count > 0){
			$flag = $flag + 1;
		}
	}
	if($flag >= 2){
		print OUT "$t[0]\t$start\t$end\t$t[3]\t$c{'plus_plus'}\t$c{'plus_minus'}\t$c{'minus_plus'}\t$c{'minus_minus'}\n";
	}
	#print Dumper(\%c);
	#print OUT "$t[0]\t$start\t$end\t$t[3]\t$c{'plus_plus'}\t$c{'plus_minus'}\t$c{'minus_plus'}\t$c{'minus_minus'}\n";
}
close IN;
close OUT;
############### merged end ##################


############################################
# filtering blacklist
`bedtools intersect -a $prefix.all.sites.merged.confirmed.raw -b $blacklist_f -v > $prefix.all.sites.merged.confirmed`;

# Filtering sites overlap with Control sample. This part is not well designed, currently.
# extend the sequence
# 2020/2/27 zhoujj2013@gmail.com
`bedtools slop -i $prefix.all.sites.merged.confirmed -g $chrom_size -b 50 > $prefix.all.sites.merged.confirmed.ext.raw.bed`;
`bedtools intersect -a $prefix.all.sites.merged.confirmed.ext.raw.bed -b $bg_f -v > $prefix.all.sites.merged.confirmed.ext.bed`;

# get fasta sequence for the targeted sites
`bedtools getfasta -fi $ref -bed $prefix.all.sites.merged.confirmed.ext.bed > $prefix.all.sites.merged.confirmed.ext.fa`;
open OUT,">","$prefix.all.sites.merged.confirmed.ext.rev.fa" || die $!;
open IN,"$prefix.all.sites.merged.confirmed.ext.fa" || die $!;
$/ = ">";<IN>;$/ = "\n";
while(<IN>){
	my $id = $1 if(/(\S+)/);
	$/ = ">";
	my $seq = <IN>;
	chomp($seq);
	$seq =~ s/\n//g;
	$/ = "\n";

	$seq = reverse($seq);
	$seq = uc($seq);
	$seq =~ tr/ATCG/TAGC/;
	print OUT ">$id\n";
	print OUT "$seq\n";
}
close IN;
close OUT;

################################################
# keep origin sgRNA for better visualization
open IN,"$grna_f" || die $!;
#open OUT,">","./grna_rev_com.fa" || die $!;
my $grna_length = 0;
my $grna;
while(<IN>){
	chomp;
	next if(/^>/);
	$grna_length = length($_);
	$grna = $_;
	#$grna = reverse($grna);
	#$grna =~ tr/TCGA/AGCT/;
	#print OUT ">grna\n";
	#print OUT "$grna\n";
}
close IN;
#close OUT;
#$grna_f = "./grna_rev_com.fa";

# index
`perl $Bin/replace_id.pl $prefix.all.sites.merged.confirmed.ext.fa > $prefix.all.sites.merged.confirmed.ext.renamed.fa 2>$prefix.all.sites.merged.confirmed.ext.renamed.index`;
`perl $Bin/replace_id.pl $prefix.all.sites.merged.confirmed.ext.rev.fa > $prefix.all.sites.merged.confirmed.ext.rev.renamed.fa 2>$prefix.all.sites.merged.confirmed.ext.rev.renamed.index`;

# align
`/usr/bin/water -asequence $grna_f -bsequence  $prefix.all.sites.merged.confirmed.ext.renamed.fa -gapopen 100 -gapextend 0.5 -outfile $prefix.align.water`;
`/usr/bin/water -asequence $grna_f -bsequence  $prefix.all.sites.merged.confirmed.ext.rev.renamed.fa -gapopen 100 -gapextend 0.5 -outfile $prefix.align.rev.water`;
########### alignment end ########################

##################################################
## process alignments with gaps
`/usr/bin/water -asequence $grna_f -bsequence  $prefix.all.sites.merged.confirmed.ext.renamed.fa -gapopen 10 -gapextend 0.5 -outfile $prefix.align.water.gap`;
`/usr/bin/water -asequence $grna_f -bsequence  $prefix.all.sites.merged.confirmed.ext.rev.renamed.fa -gapopen 10 -gapextend 0.5 -outfile $prefix.align.rev.water.gap`;

# parse gap information
`perl $Bin/parsing_water_for_gap_mapping.v2.pl $prefix.align.water.gap $prefix.align.rev.water.gap $prefix.all.sites.merged.confirmed.ext.renamed.fa $prefix.all.sites.merged.confirmed.ext.rev.renamed.fa $prefix.all.sites.merged.confirmed.ext.renamed.index $prefix.all.sites.merged.confirmed.ext.rev.renamed.index $ref $grna $max_mismatch $max_gap > ./$prefix.parsing_water_for_visualization.offtarget.gap.raw.bed`;
`perl $Bin/get_mapped_sites4gaps.pl $prefix.parsing_water_for_visualization.offtarget.gap.raw.bed $prefix.all.sites.merged.confirmed.ext.fa $prefix.all.sites.merged.confirmed.ext.rev.fa $prefix.all.sites.merged.confirmed.ext.bed > $prefix.parsing_water_for_visualization.offtarget.gap.bed`;
########### end ###################################

###################################################
## retrive alignment information and recover the coordinates
#`perl /home/zhoujj/project/02guide_seq/bin/parsing_water.pl $prefix.align.water $prefix.all.sites.merged.confirmed.ext.renamed.index > $prefix.align.water.result`;
#`perl /home/zhoujj/project/02guide_seq/bin/parsing_water.pl $prefix.align.rev.water $prefix.all.sites.merged.confirmed.ext.rev.renamed.index > $prefix.align.water.rev.result`;
#
### recover the coordinates
#`perl /home/zhoujj/project/02guide_seq/bin/recover_coor.pl $prefix.align.water.result $prefix.all.sites.merged.confirmed.ext.bed $prefix.all.sites.merged.confirmed.ext.renamed.index > $prefix.clustered.region.ext.align.bed`;
#`perl /home/zhoujj/project/02guide_seq/bin/recover_coor.pl $prefix.align.water.result $prefix.all.sites.merged.confirmed.ext.bed $prefix.all.sites.merged.confirmed.ext.rev.renamed.index > $prefix.clustered.region.ext.rev.align.bed`;
############## the old version ####################

###################################################
# process alignments without gaps 
`perl $Bin/parsing_water_for_visualization.pl $prefix.align.water $prefix.align.rev.water $prefix.all.sites.merged.confirmed.ext.renamed.fa $prefix.all.sites.merged.confirmed.ext.rev.renamed.fa $prefix.all.sites.merged.confirmed.ext.renamed.index $prefix.all.sites.merged.confirmed.ext.rev.renamed.index $grna $grna_length $max_mismatch > $prefix.parsing_water_for_visualization.txt`;

`perl $Bin/prepare_for_visualization.v2.pl $prefix.all.sites.merged.confirmed.ext.renamed.index $prefix.all.sites.merged.confirmed.ext.bed $prefix.parsing_water_for_visualization.txt > $prefix.parsing_water_for_visualization.offtarget.raw`;

# merge identical target site
`perl $Bin/merge_identical_sites.pl $prefix.parsing_water_for_visualization.offtarget.raw > $prefix.parsing_water_for_visualization.offtarget`;

`sed '1d' $prefix.parsing_water_for_visualization.offtarget | cut -f 1-3,4,12,22,33,34 | bedtools sort -i - | sort -k5nr > $prefix.parsing_water_for_visualization.offtarget.bed`;
################# end ############################

##################################################
# combine the nongap alignment and gap alignment
##################################################
# for record the detail information
`perl $Bin/combine_nongaps_gaps.pl $prefix.parsing_water_for_visualization.offtarget.bed $prefix.parsing_water_for_visualization.offtarget.gap.bed > $prefix.parsing_water_for_visualization.offtarget.combined.bed`;

# for final visualization
#print STDERR "perl $Bin/generate_gap_visualization.pl ./$prefix.parsing_water_for_visualization.offtarget.gap.raw.bed ./$prefix.parsing_water_for_visualization.offtarget.gap.bed ./$prefix.parsing_water_for_visualization.offtarget.bed > ./$prefix.parsing_water_for_visualization.offtarget.combined.forVis.bed\n";

`perl $Bin/generate_gap_visualization.pl ./$prefix.parsing_water_for_visualization.offtarget.gap.raw.bed ./$prefix.parsing_water_for_visualization.offtarget.gap.bed ./$prefix.parsing_water_for_visualization.offtarget.bed > ./$prefix.parsing_water_for_visualization.offtarget.combined.forVis.bed`;
################### end ##########################


################ visualization ###################
mkdir "./$prefix.drawTargetsite" unless(-d "./$prefix.drawTargetsite");

#print "python /home/zhoujj/software/guideseq/guideseq/guideseq.py visualize --infile $prefix.parsing_water_for_visualization.offtarget --outfolder $prefix.drawTargetsite/ > $prefix.drawTargetsite/$prefix.log 2> $prefix.rawTargetsite/$prefix.err\n";

#print "python $Bin/visualization.py --identified_file $prefix.parsing_water_for_visualization.offtarget.combined.forVis.bed --outfile $prefix.drawTargetsite/$prefix\_offtargets --title $prefix --PAM $pam\n";

`python $Bin/visualization.py --identified_file $prefix.parsing_water_for_visualization.offtarget.combined.forVis.bed --outfile $prefix.drawTargetsite/$prefix\_offtargets --title $prefix --PAM $pam`;
`rsvg-convert -f pdf -o $prefix.drawTargetsite/$prefix\_offtargets.pdf $prefix.drawTargetsite/$prefix\_offtargets.svg`;
##################### end #######################
