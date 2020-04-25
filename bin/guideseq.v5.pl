#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);

&usage if @ARGV<1;

sub usage {
        my $usage = << "USAGE";

        This script is designed for TAG-seq analysis (not including DNA-seq based variant calling).
        Author: zhoujj2013\@gmail.com 
        Usage: $0 <config.txt>

USAGE
print "$usage";
exit(1);
};

my $conf=shift;
my %conf;
&load_conf($conf, \%conf);

my $all='all: ';
my $mk;

my $sample_f = abs_path($conf{SAMPLE});
my $out = abs_path($conf{OUTDIR});
my $thread = $conf{THREAD};

#### write you things ###
#decide to add function in another file
#my %compare;
open IN,"$sample_f" || die $!;
while(<IN>){
	chomp;
	next if(/^$/);
	next if(/^#/);
	my @t = split /\t/;
	my $id = shift @t;
	##################################################
	# decide to add function in another file
	##################################################
	#my $control_case = pop(@t); 
	#
	my $tag = pop(@t);
	my $seq_strand = pop(@t);
	my $tag_rev_comp = reverse($tag);
	$tag_rev_comp =~ tr/TCGA/AGCT/;
	#print "$tag_rev_comp\n";
	##################################################
        # decide to add function in another file
        ##################################################
	#push @{$compare{$control_case}},$id;

	# determine PE
	my $type = "PE";
	my ($r1,$r2);
	if(scalar(@t) == 1){
		$type ="SE";
		$r1 = abs_path($t[0]);
		$r2 = $r1;
	}elsif(scalar(@t) == 2){
		$type = "PE";
		$r1 = abs_path($t[0]);
		$r2 = abs_path($t[1]);
	}
		
	## QC
	make_path "$out/$id/00datafilter";
	my $encode_type = `perl $conf{BIN}/phredDetector.pl $r1`;
	my $filter_cutoff = 5;
	my $phred = "";
    	if($encode_type == 33){
        	#$filter_cutoff = 38;
        	$phred = " ";
    	}elsif($encode_type == 64){
        	#$filter_cutoff = 69;
        	$phred = " --phred64 ";
    	}

	$mk .= "00fastqc.origin.finished: $r1 $r2\n";
	$mk .= "\t$conf{FASTQC} --threads $conf{THREAD} -f fastq -o $out/$id/00datafilter/ $r1 $r2 > /dev/null 2>/dev/null && touch 00fastqc.origin.finished\n";
	$all .= "00fastqc.origin.finished ";

	# remove ODN
	$mk .= "00rmODN.finished: $r1 $r2\n";
	$mk .= "\tperl $conf{BIN}/remove_ODN.pl $r1 $r2 $tag_rev_comp $out/$id/00datafilter $id > $out/$id/00datafilter/$id.rmODN.log 2> $out/$id/00datafilter/$id.rmODN.err && touch 00rmODN.finished\n";
	#$mk .= "\tperl $conf{BIN}/remove_ODN.pl $out/$id/00datafilter/$id.Trim.R1.fq $out/$id/00datafilter/$id.Trim.R2.fq $tag_rev_comp > $out/$id/00datafilter/$id.Trim.rmODN.R2.fq 2> $out/$id/00datafilter/$id.Trim.rmODN.R1.fq && touch 00rmODN.finished\n";
	$all .= "00rmODN.finished ";
	
	$r1 = "$out/$id/00datafilter/$id.rmODN.R1.fq";
	$r2 = "$out/$id/00datafilter/$id.rmODN.R2.fq";

	$mk .= "00trim.finished: 00rmODN.finished\n";	
	$mk .= "\t$conf{AdapterRemoval} --qualitybase $encode_type  --file1 $r1 --file2 $r2 --settings  $out/$id/00datafilter/$id.adapterRemoval.setting --output1 $out/$id/00datafilter/$id.Trim.R1.fq --output2 $out/$id/00datafilter/$id.Trim.R2.fq --singleton $out/$id/00datafilter/$id.Single.fq --discarded $out/$id/00datafilter/$id.Discarded.fq --minlength $conf{MINLEN} --trimns --trimqualities --minquality $filter_cutoff --trimwindows 15 --threads $conf{THREAD} --adapter-list $conf{ADAPTER} > $out/$id/00datafilter/$id.trim.log 2>$out/$id/00datafilter/$id.trim.err && touch 00trim.finished\n";
	#$mk .= "\t$Bin/Ktrim.v2.3 -a $r1 -b $r2 -o $out/$id/00datafilter/$id.Trim -5 AAAAAAAAAAAAAAAAAAAAAAAAAAA -3 AAAAAAAAAAAAAAAAAAAAAAAAAAA > $out/$id/00datafilter/$id.trim.log 2>$out/$id/00datafilter/$id.trim.err && touch 00trim.finished\n";
	$all .= "00trim.finished ";

	# remove UMIs
	$mk .= "00rmumis.finished: 00trim.finished\n";
	$mk .= "\t$conf{umi_tools} extract -I $out/$id/00datafilter/$id.Trim.R1.fq --bc-pattern=NNNNNNNN --read2-in=$out/$id/00datafilter/$id.Trim.R2.fq --stdout=$out/$id/00datafilter/$id.Trim.R1.umis.fq --read2-out=$out/$id/00datafilter/$id.Trim.R2.umis.fq > $out/$id/00datafilter/$id.Trim.umis.log 2>$out/$id/00datafilter/$id.Trim.umis.err && touch 00rmumis.finished\n";
	$all .= "00rmumis.finished ";

	$mk .= "00fastqc.finished: 00rmumis.finished\n";
	$mk .= "\t$conf{FASTQC} --threads $conf{THREAD} -f fastq -o $out/$id/00datafilter/ $out/$id/00datafilter/$id.Trim.R1.fq $out/$id/00datafilter/$id.Trim.R2.fq > /dev/null 2>/dev/null && touch 00fastqc.finished\n";
	$all .= "00fastqc.finished ";
	
	#if($seq_strand eq "forward"){
	#	$mk .= "00strip.finished: 00trim.finished\n";
	#	$mk .= "\tperl $conf{BIN}/forward.pl $out/$id/00datafilter/$id.Trim.R2.fq > $out/$id/00datafilter/$id.Trim.strip.fq 2>/dev/null && perl $conf{BIN}/count_fastq_len.pl $out/$id/00datafilter/$id.Trim.strip.fq > $out/$id/00datafilter/$id.insert_size && Rscript $conf{BIN}/count_fastq_len_plot.r $out/$id/00datafilter/$id.insert_size && touch 00strip.finished\n";
	#	$all .= "00strip.finished ";
	#}elsif($seq_strand eq "reverse"){
	#	$mk .= "00strip.finished: 00trim.finished\n";
	#	$mk .= "\tperl $conf{BIN}/reverse.pl $out/$id/00datafilter/$id.Trim.R2.fq > $out/$id/00datafilter/$id.Trim.strip.fq 2>/dev/null && perl $conf{BIN}/count_fastq_len.pl $out/$id/00datafilter/$id.Trim.strip.fq > $out/$id/00datafilter/$id.insert_size && Rscript $conf{BIN}/count_fastq_len_plot.r $out/$id/00datafilter/$id.insert_size && touch 00strip.finished\n";
	#	$all .= "00strip.finished ";
	#}
		
	## alignment
	#my $single = "$out/$id/00datafilter/$id.Single.fq";
	#my $discarded = "$out/$id/00datafilter/$id.Discarded.fq";
	my $clear_r1 = "$out/$id/00datafilter/$id.Trim.R1.umis.fq";
	my $clear_r2 = "$out/$id/00datafilter/$id.Trim.R2.umis.fq";
	
	make_path "$out/$id/01alignment";
	$mk .= "01align.finished: 00fastqc.finished\n";
	# --outSAMtype BAM Unsorted ###SortedByCoordinate
	# add --alignIntronMax 1 to avoid large intron mapping
	# take care of --alignIntronMax 50 --outFilterScoreMinOverLread 0.5
	# this is a risk parameters
	#
	$mk .= "\t$conf{STAR} --genomeDir $conf{INDEX} --runThreadN $thread --readFilesIn $clear_r1 $clear_r2 --outFileNamePrefix $out/$id/01alignment/$id. --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --alignIntronMax 50 --outFilterScoreMinOverLread 0.5 1>$out/$id/01alignment/$id.star.log 2>$out/$id/01alignment/$id.star.err && $conf{SAMTOOLS} index $out/$id/01alignment/$id.Aligned.sortedByCoord.out.bam && touch 01align.finished\n";
	$all .= "01align.finished ";

	$mk .= "01rmdup.finished: 01align.finished\n";
	# rational: https://cgatoxford.wordpress.com/2015/08/14/unique-molecular-identifiers-the-problem-the-solution-and-the-proof/
	# https://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract
	$mk .= "\t$conf{umi_tools} dedup -I $out/$id/01alignment/$id.Aligned.sortedByCoord.out.bam --paired -S $out/$id/01alignment/$id.Aligned.sortedByCoord.out.dedup.bam --output-stats=$out/$id/01alignment/$id.dedup > $out/$id/01alignment/$id.dedup.log 2>$out/$id/01alignment/$id.dedup.err && $conf{SAMTOOLS} flagstat $out/$id/01alignment/$id.Aligned.sortedByCoord.out.dedup.bam > $out/$id/01alignment/$id.Aligned.sortedByCoord.out.dedup.bam.flagstat && touch 01rmdup.finished\n";
	#$mk .= "\tjava -jar $conf{PICARD} MarkDuplicates I=$out/$id/01alignment/$id.Aligned.sortedByCoord.out.bam O=$out/$id/01alignment/$id.sorted.rmdup.bam M=$out/$id/01alignment/$id.rmdup.log REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR=./ >$out/$id/01alignment/rmdup.log 2>$out/$id/01alignment/rmdup.err && $conf{SAMTOOLS} index $out/$id/01alignment/$id.sorted.rmdup.bam && touch 01rmdup.finished\n";
	$all .= "01rmdup.finished ";

	make_path "$out/$id/02potentialTargets";
	$mk .= "02potentialTargets.finished: 01rmdup.finished\n";
	# ref: https://broadinstitute.github.io/picard/explain-flags.html
	# add -q 30 to filter low quality mapping in repeat region
	# remove perl $conf{BIN}/filter8odn.v4.pl $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.sam $tag $seq_strand > $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.filtered.sam
	# NOTE: import setting for detect peaks, we use 10bp windows for calculating coverage, the minimue distance for detect peaks, this method can avoid replicately calculate the coverage. 2019/12/11
	$mk .= "\t$conf{SAMTOOLS} view -f 128 -q 30 $out/$id/01alignment/$id.Aligned.sortedByCoord.out.dedup.bam > $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.filtered.sam && $conf{SAMTOOLS} view -H $out/$id/01alignment/$id.Aligned.sortedByCoord.out.bam | cat - $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.filtered.sam | $conf{SAMTOOLS} view -b - > $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.filtered.bam && $conf{SAMTOOLS} flagstat $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.filtered.bam > $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.filtered.bam.flagstat && $conf{BEDTOOLS} bamtobed -i $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.filtered.bam > $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.filtered.bed && perl $conf{BIN}/extractLoci.pl $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.filtered.bed > $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.filtered.loci.bed && $conf{BEDTOOLS} sort -i $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.filtered.loci.bed > $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.filtered.loci.sorted.bed && $conf{BEDTOOLS} merge -d 10 -c 4 -o count -i $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.filtered.loci.sorted.bed > $out/$id/02potentialTargets/$id.clustered.region.bed && $conf{BEDTOOLS} slop -i $out/$id/02potentialTargets/$id.clustered.region.bed -g $conf{CHROMSIZE} -b 10 > $out/$id/02potentialTargets/$id.clustered.region.extend10bp.bed && perl $conf{BIN}/replace_id_bed.pl $out/$id/02potentialTargets/$id.clustered.region.extend10bp.bed $id > $out/$id/02potentialTargets/$id.clustered.region.rpid.bed && $conf{BEDTOOLS} makewindows -b $out/$id/02potentialTargets/$id.clustered.region.rpid.bed -w 10 -s 1 -i srcwinnum > $out/$id/02potentialTargets/$id.clustered.region.rpid.wins.bed && perl $conf{BIN}/seperate_plus_minus.pl $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.filtered.loci.sorted.bed > $out/$id/02potentialTargets/$id.plus.bed 2>$out/$id/02potentialTargets/$id.minus.bed && $conf{bedops}/bedmap --echo --count --delim \'\\t\' $out/$id/02potentialTargets/$id.clustered.region.rpid.wins.bed $out/$id/02potentialTargets/$id.plus.bed > $out/$id/02potentialTargets/$id.plus.count.raw && perl $conf{BIN}/add_strand_flag.pl $out/$id/02potentialTargets/$id.plus.count.raw plus > $out/$id/02potentialTargets/$id.plus.count && $conf{bedops}/bedmap --echo --count --delim \'\\t\' $out/$id/02potentialTargets/$id.clustered.region.rpid.wins.bed $out/$id/02potentialTargets/$id.minus.bed > $out/$id/02potentialTargets/$id.minus.count.raw && perl $conf{BIN}/add_strand_flag.pl $out/$id/02potentialTargets/$id.minus.count.raw minus > $out/$id/02potentialTargets/$id.minus.count && python $conf{BIN}/detect.peaks.py $out/$id/02potentialTargets/$id.plus.count $out/$id/02potentialTargets/$id.minus.count 5 10 >$out/$id/02potentialTargets/$id.plus.proximal 2>$out/$id/02potentialTargets/$id.minus.proximal && $conf{BEDTOOLS} sort -i $out/$id/02potentialTargets/$id.plus.proximal > $out/$id/02potentialTargets/$id.plus.proximal.sorted.raw && $conf{BEDTOOLS} sort -i $out/$id/02potentialTargets/$id.minus.proximal > $out/$id/02potentialTargets/$id.minus.proximal.sorted.raw && perl $conf{BIN}/strand_nomalization.pl $out/$id/02potentialTargets/$id.minus.proximal.sorted.raw minus > $out/$id/02potentialTargets/$id.minus.proximal.sorted && perl $conf{BIN}/strand_nomalization.pl $out/$id/02potentialTargets/$id.plus.proximal.sorted.raw plus > $out/$id/02potentialTargets/$id.plus.proximal.sorted && touch 02potentialTargets.finished\n";
	########### the result is stable ##########
	# 5	102
	# 10	101
	# 20	99
	# 30	99
	# 40	99
	# 50	98
	# 100	97
	# #########################################
	$all .= "02potentialTargets.finished ";

	# findTarget script has obsoleted. 2020/3/3
	#make_path "$out/$id/03findTargets";
	#$mk .= "03findTargets.finished: 02potentialTargets.finished\n";
	#$mk .= "\t$conf{BEDTOOLS} slop -i $out/$id/02potentialTargets/$id.clustered.region.bed -g $conf{hg19}/hg19.chrom.sizes -b 25 > $out/$id/03findTargets/$id.clustered.region.ext25.bed && $conf{BEDTOOLS} getfasta -fi $conf{hg19}/hg19.fa -bed $out/$id/03findTargets/$id.clustered.region.ext25.bed > $out/$id/03findTargets/$id.clustered.region.ext25.fa && perl $conf{BIN}/replace_id.pl $out/$id/03findTargets/$id.clustered.region.ext25.fa > $out/$id/03findTargets/$id.clustered.region.ext25.renamed.fa 2>$out/$id/03findTargets/$id.clustered.region.ext25.renamed.index && $conf{water} -asequence $conf{GRNA} -bsequence $out/$id/03findTargets/$id.clustered.region.ext25.renamed.fa -gapopen 10 -gapextend 0.5 -outfile $out/$id/03findTargets/$id.align.water && perl $conf{BIN}/parsing_water.pl $out/$id/03findTargets/$id.align.water > $out/$id/03findTargets/$id.align.water.result && perl $conf{BIN}/recover_coor.pl $out/$id/03findTargets/$id.align.water.result $out/$id/03findTargets/$id.clustered.region.ext25.bed $out/$id/03findTargets/$id.clustered.region.ext25.renamed.index > $out/$id/03findTargets/$id.clustered.region.ext25.align.bed && touch 03findTargets.finished\n";
	#$all .= "03findTargets.finished ";

	make_path "$out/$id/03visual";
	$mk .= "03visual.finished: 02potentialTargets.finished\n";
	$mk .= "\tperl $conf{BIN}/loci2bdg.pl $out/$id/02potentialTargets/$id.Aligned.sortedByCoord.out.dedup.single.filtered.loci.bed > $out/$id/03visual/$id.plus.bdg 2> $out/$id/03visual/$id.minus.bdg && $conf{BEDTOOLS} sort -i $out/$id/03visual/$id.minus.bdg > $out/$id/03visual/$id.minus.sorted.bdg && $conf{BEDTOOLS} sort -i $out/$id/03visual/$id.plus.bdg > $out/$id/03visual/$id.plus.sorted.bdg && $conf{BIN}/bedGraphToBigWig $out/$id/03visual/$id.minus.sorted.bdg $conf{CHROMSIZE} $out/$id/03visual/$id.minus.bw > $out/$id/03visual/$id.minus.bw.log 2> $out/$id/03visual/$id.minus.bw.err && $conf{BIN}/bedGraphToBigWig $out/$id/03visual/$id.plus.sorted.bdg $conf{CHROMSIZE} $out/$id/03visual/$id.plus.bw > $out/$id/03visual/$id.plus.bw.log 2> $out/$id/03visual/$id.plus.bw.err && touch 03visual.finished\n";
	$all .= "03visual.finished ";

	## 04checkDSBs script has obsoleted, because the cell type speicific of DSBs. 2020/3/3
	#make_path "$out/$id/04checkDSBs";
	#$mk .= "04checkDSBs.finished: 03visual.finished\n";
	#$mk .= "\t$conf{BEDTOOLS} intersect -a $out/$id/03findTargets/$id.clustered.region.ext25.align.bed -b $conf{BIN}/data/Miotic_DSBs.bed > $out/$id/04checkDSBs/$id.ol.Miotic_DSBs.bed && $conf{BEDTOOLS} intersect -a $out/$id/03findTargets/$id.clustered.region.ext25.align.bed -b $conf{BIN}/data/spontaneous_DSBs.bed > $out/$id/04checkDSBs/$id.ol.spontaneous_DSBs.bed && touch 04checkDSBs.finished\n";
	#$all .= "04checkDSBs.finished ";

	#make_path "$out/$id/06checkOverlap";
	#$mk .= "06checkOverlap.finished: 04checkDSBs.finished\n";
	#$mk .= "\t$conf{BEDTOOLS} intersect -a $out/$id/03findTargets/$id.clustered.region.ext25.align.bed -b $out/$id/04visual/$id.plus.sorted.bdg -wo > $out/$id/06checkOverlap/$id.potential.target.vs.plus.intersect && $conf{BEDTOOLS} intersect -a $out/$id/03findTargets/$id.clustered.region.ext25.align.bed -b $out/$id/04visual/$id.minus.sorted.bdg -wo > $out/$id/06checkOverlap/$id.potential.target.vs.minus.intersect && touch 06checkOverlap.finished\n";
	#$all .= "06checkOverlap.finished ";

	#make_path "$out/$id/07drawTargetsite";
	#$mk .= "07drawTargetsite.finished: 06checkOverlap.finished\n";
	#$mk .= "\tperl $conf{BIN}/parsing_water_for_visualization.pl $out/$id/03findTargets/$id.align.water $conf{GRNA_LENGTH} > $out/$id/07drawTargetsite/$id.parsing_water_for_visualization.txt && perl $conf{BIN}/prepare_for_visualization.pl $out/$id/03findTargets/$id.clustered.region.ext25.renamed.index $out/$id/03findTargets/$id.clustered.region.ext25.bed $out/$id/07drawTargetsite/$id.parsing_water_for_visualization.txt > $out/$id/07drawTargetsite/$id.parsing_water_for_visualization.offtarget && python $conf{guideseq}/guideseq/guideseq.py visualize --infile $out/$id/07drawTargetsite/$id.parsing_water_for_visualization.offtarget --outfolder $out/$id/07drawTargetsite/ > $out/$id/07drawTargetsite/$id.log 2> $out/$id/07drawTargetsite/$id.err && touch 07drawTargetsite.finished\n";
	#$all .= "07drawTargetsite.finished ";
	
	make_path abs_path($conf{OUTDIR});
	open OUT, ">$out/$id/makefile";
	print OUT $all, "\n";
	print OUT $mk, "\n";
	close OUT;
	$all = "all: ";
	$mk = "";
}
close IN;

# write the comparing part



#########################
sub load_conf
{
    my $conf_file=shift;
    my $conf_hash=shift; #hash ref
    open CONF, $conf_file || die "$!";
    while(<CONF>)
    {
        chomp;
        next unless $_ =~ /\S+/;
        next if $_ =~ /^#/;
        warn "$_\n";
        my @F = split"\t", $_;  #key->value
        $conf_hash->{$F[0]} = $F[1];
    }
    close CONF;
}
