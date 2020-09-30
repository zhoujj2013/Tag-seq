#!/usr/bin/perl -w

use strict;

&usage if @ARGV<7;

sub usage {
        my $usage = << "USAGE";

        Detect gap open alignments.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 config.cfg

USAGE
print "$usage";
exit(1);
};

#############################################################
# input files
my $fw = shift; # waterman alignment result, *.water.gap
my $rv = shift; # waterman alignment result, *.rev.water.gap

my $fw_fa = shift; # input renamed fa, forward strand
my $rv_fa = shift; # input renamed fa, reverse strand
my $fw_renamed_index_f = shift;
my $rv_renamed_index_f = shift;

my $hg19_f = shift;
my $grna = shift; # grna fasta file
my $grna_len = length($grna);

# change the mismatch and gaps
#my $pre_mismatch_and_gap = shift;
my $gap_number_cutoff = shift;
my $gap_mismatch_cutoff = shift;

#############################################################

#############################################################
# read in the index file, so that we can trace coordinates.
# start
my %index_fw;
my %index_rv;

open IN,"$fw_renamed_index_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	$index_fw{$t[0]} = $t[1];
}
close IN;

open IN,"$rv_renamed_index_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	$index_rv{$t[0]} = $t[1];
}
close IN;
##############################################################

##############################################################
# Read in sequence.
# Start
my %fwseq;
open IN,"$fw_fa" || die $!;
$/ = ">";<IN>;$/ = "\n";
while(<IN>){
	my $id = $1 if(/(\S+)/);
	$/ = ">";
	my $seq = <IN>;
	chomp($seq);
	$seq =~ s/\n//g;
	$/ = "\n";
	$seq = uc($seq);
	$fwseq{$id} = $seq;
}
close IN;

my %rvseq;
open IN,"$rv_fa" || die $!;
$/ = ">";<IN>;$/ = "\n";
while(<IN>){
	my $id = $1 if(/(\S+)/);
	$/ = ">";
	my $seq = <IN>;
	chomp($seq);
	$seq =~ s/\n//g;
	$/ = "\n";
	$seq = uc($seq);
	$rvseq{$id} = $seq;
}
close IN;

my %hg19;
open IN,"$hg19_f" || die $!;
$/ = ">";<IN>;$/ = "\n";
while(<IN>){
	chomp;
	my $id = $1 if(/(\S+)/);
	$/ = ">";
	my $seq = <IN>;
	chomp($seq);
	$seq =~ s/\n//g;
	$/ = "\n";
	$seq = uc($seq);
	$hg19{$id} = $seq;
}
close IN;

# end
#################################################################

#print "## Forward strand\n";
#print "Ref\tOT-REF $grna\n";
# Process water alignment, 
open IN,"$fw" || die $!;
$/ = "#======================================="; <IN>; $/ = "\n";
while(<IN>){
	$/ = "#=======================================";
	my $con = <IN>;
	chomp($con);
	my $con2 = <IN>;
	my $align = $con2;
	chomp($align);
	$/ = "\n";

	# recover coordinates
	#######################
	##=======================================
	#
	#grna               1 CCCCATCATGTAGGTT     16
	#                     ||.||..|.|.||.||
	#1                  7 cctcacaaaggagttt     22
	#
	my @ali = split /\n/,$align;
	my ($s,$e) = ($2,$3) if($ali[4] =~ /(\d+)\s+(\d+)\s+\S+\s+(\d+)/);
	#print "$s,$e\n";

	my $seq1 = $1 if($con =~ /# 1: (\S+)/);
	my $seq2 = $1 if($con =~ /# 2: (\S+)/);
	
	my $length = $1 if($con =~ /# Length: (\S+)/);
	my $identity = $1 if($con =~ /# Identity:.*\((\S+)\)/);
	my $similarity = $1 if($con =~ /# Similarity:.*\((\S+)\)/);

	# get gap information
	my $gap = 0; # gap percentage
	my $gap2 = 0; # gap number (absoluated number)
	if($con =~ /# Gaps:\s+(\d+)\/.*\((.*)\)/){
		$gap2 = $1;
		$gap = $2;
		$gap =~ s/%//g;
		$gap =~ s/ //g;
		$gap =~ s/\s+//g;

	}
	#print "$gap2\n";
	# skin alignment without gaps
	next if($gap2 < 1 || $gap2 >= 2);

	my ($gs,$gseq,$ge) = ($1,$2,$3) if($ali[2] =~ /grna\s+(\d+)\s+(\S+)\s+(\d+)/);
	my ($rid,$rs,$rseq,$re) = ($1,$2,$3,$4) if($ali[4] =~ /(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/);
	my ($ali_str) = ($1) if($ali[3] =~ /\s+(.*)/); # align string: ||.||..|.|.||.||

	# recover the coordinates of reference, 1 base coordinate system.
	my @coor = ($1,$2,$3) if($index_fw{$rid} =~ /(\w+):(\d+)-(\d+)/);
	my $ors = $coor[1] + $rs;
	my $ore = $ors + $re - $rs;
	my @new_coor = ($coor[0],$ors,$ore);

	#print "##############aabb################\n";
	#my $tseq = getSeq(\@new_coor,\%hg19);
	#print "$align\n";
	#print "$gseq\n";
	#print "$ali_str\n";
	#print "$rseq\n";
	#print "$coor[0]\t$ors\t$ore\n";
	#print "$gseq\n";
	#print "$rseq\n";
	#print "$tseq\n";
	#print "##############################\n";

	# get end to end alignment, because water is local alignment software.
	my $g_left_s = 1;
	my $g_left_e = $gs-1;
	my $g_right_s = $ge + 1;
	my $g_right_e = $grna_len;
	
	my $gseq_l = "";
	my $gseq_r = "";
	if($g_left_e < 1){
		$gseq_l = "";
	}else{
		$gseq_l = substr($grna, $g_left_s - 1, $g_left_e - $g_left_s + 1);
	}

	if($g_right_e < $g_right_s){
		$gseq_r = "";
	}else{
		$gseq_r = substr($grna, $g_right_s - 1, $g_right_e - $g_right_s + 1);
	}

	my $r_left_s = 0;
	my $r_left_e = 0;
	my $r_right_s = 0;
	my $r_right_e = 0;

	my $tseq_l = "";
	my $tseq_r = "";
	
	if($g_left_e < 1){
		$tseq_l = "";
	}else{
		$r_left_e = $ors-1;
		$r_left_s = $ors-($g_left_e-$g_left_s+1);
		my @r_left = ($coor[0],$r_left_s,$r_left_e);
		$tseq_l = getSeq(\@r_left,\%hg19);
	}

	if($g_right_e < $g_right_s){
		$tseq_r = "";
	}else{
		$r_right_s = $ore + 1;
		$r_right_e = $ore + ($g_right_e - $g_right_s + 1);
		my @r_right = ($coor[0],$r_right_s,$r_right_e);
		$tseq_r = getSeq(\@r_right,\%hg19);
	}
	#######################################################
	# Get the right coordinate, 2020.0.26
	#######################################################
	
	#print "#########$ali_str#########\n";
	#print "$g_left_s\t$g_left_e\t$g_right_s\t$g_right_e\n";
	#print "$r_left_s\t$r_left_e\t$r_right_s\t$r_right_e\n";
	#print "$grna\n";

	my $realign = "$gseq_l"."$gseq"."$gseq_r";
	my $bugle = "$tseq_l"."$rseq"."$tseq_r";
	$realign =~ tr/atgc/ATGC/;
	$bugle =~ tr/atgc/ATGC/;

	#print "$realign\n";
	#print "$bugle\n";
	#print "#############cccdddd#################\n";

	my $chr_id = $coor[0];
	my $r_start = $coor[1];
	my $r_end = $coor[2];
	my $gap_number = 0;
	my $mis = 0;

	($mis, $gap_number) = align($realign,$bugle);
	#if($gap_number <= $gap_number_cutoff && $mis + $gap_number <= $pre_mismatch_and_gap){
	if($gap_number <= $gap_number_cutoff && $mis <= $gap_mismatch_cutoff){
		print "$chr_id\t$r_start\t$r_end\t$gap_number\t$mis\t$grna\t$realign\t$ali_str\t$bugle\t$rid\t+\n";
	}
}
close IN;

#print "# Reverse strand\n";
## process water alignment (reverse part)
open IN,"$rv" || die $!;
$/ = "#======================================="; <IN>; $/ = "\n";
while(<IN>){
	$/ = "#=======================================";
	my $con = <IN>;
	chomp($con);
	my $con2 = <IN>;
	my $align = $con2;
	chomp($align);
	$/ = "\n";

	# recover coordinates
	my @ali = split /\n/,$align;
	my ($s,$e) = ($2,$3) if($ali[4] =~ /(\d+)\s+(\d+)\s+\S+\s+(\d+)/);

	my $seq1 = $1 if($con =~ /# 1: (\S+)/);
	my $seq2 = $1 if($con =~ /# 2: (\S+)/);

	my ($chr,$start,$end) = ($1, $2, $3) if($index_rv{$seq2} =~ /(.*):(\d+)-(\d+)/);
	my $new_start = $end - $e;
	my $new_end = $end - $s + 1;

	my $length = $1 if($con =~ /# Length: (\S+)/);
	my $identity = $1 if($con =~ /# Identity:.*\((\S+)\)/);
	my $similarity = $1 if($con =~ /# Similarity:.*\((\S+)\)/);
	
	# get gap information
	my $gap = 0;
	my $gap2 = 0;
	if($con =~ /# Gaps:\s+(\d+)\/.*\((.*)\)/){
		$gap2 = $1;
		$gap = $2;
		$gap =~ s/%//g;
		$gap =~ s/ //g;
		$gap =~ s/\s+//g;
	}

	# skip alignment without gaps
        next if($gap2 < 1 || $gap2 >= 2);
	
	my ($gs,$gseq,$ge) = ($1,$2,$3) if($ali[2] =~ /grna\s+(\d+)\s+(\S+)\s+(\d+)/);
	my ($rid,$rs,$rseq,$re) = ($1,$2,$3,$4) if($ali[4] =~ /(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/);
	my ($ali_str) = ($1) if($ali[3] =~ /\s+(.*)/);
	
	# recover the coordinates of reference, 1 base coordinate system.
	my @coor = ($1,$2,$3) if($index_rv{$rid} =~ /(\w+):(\d+)-(\d+)/);
	my $ore = $coor[2] - $rs + 1;
	my $ors = $ore - ($re - $rs); #$coor[1] + $rs;
	#my $ore = $coor[2] - $rs; #+ $re - $rs; #$ors + $re - $rs;
	my @new_coor = ($coor[0],$ors,$ore);

	#print "##############aabb################\n";
	#my $tseq = getSeq(\@new_coor,\%hg19);
	#$tseq = reverse_comp($tseq);=
	#print "$align\n";
	#print "$gseq\n";
	#print "$ali_str\n";
	#print "$rseq\n";
	#print "$coor[0]\t$ors\t$ore\n";
	#print "$gseq\n";
	#print "$rseq\n";
	#print "$tseq\n";
	#print "##############################\n";

	# get end to end alignment, because water is local alignment software.	
	my $g_left_s = 1;
	my $g_left_e = $gs-1;
	my $g_right_s = $ge + 1;
	my $g_right_e = $grna_len;
	
	my $gseq_l = "";
	my $gseq_r = "";
	if($g_left_e < 1){
		$gseq_l = "";
	}else{
		$gseq_l = substr($grna, $g_left_s - 1, $g_left_e - $g_left_s + 1);
	}

	if($g_right_e < $g_right_s){
		$gseq_r = "";
	}else{
		$gseq_r = substr($grna, $g_right_s - 1, $g_right_e - $g_right_s + 1);
	}

	my $r_left_s = 0;
	my $r_left_e = 0;
	my $r_right_s = 0;
	my $r_right_e = 0;

	my $tseq_l = "";
	my $tseq_r = "";
	
	if($g_left_e < 1){
		$tseq_l = "";
	}else{
		$r_left_s = $ore + 1;
		$r_left_e = $ore + ($g_left_e - $g_left_s + 1);
		my @r_left = ($coor[0],$r_left_s,$r_left_e);
		$tseq_l = getSeq(\@r_left,\%hg19);
		$tseq_l = reverse_comp($tseq_l);
	}

	if($g_right_e < $g_right_s){
		$tseq_r = "";
	}else{
		$r_right_e = $ors - 1;
		$r_right_s = $ors - ($g_right_e - $g_right_s + 1);
		my @r_right = ($coor[0],$r_right_s,$r_right_e);
		$tseq_r = getSeq(\@r_right,\%hg19);
		$tseq_r = reverse_comp($tseq_r);
	}

	#######################################################
	# Get the right coordinate, 2020.0.26
	#######################################################
	#	
	#print "#########$ali_str#########\n";
	#print "$g_left_s\t$g_left_e\t$g_right_s\t$g_right_e\n";
	#print "$r_left_s\t$r_left_e\t$r_right_s\t$r_right_e\n";
	#print "$grna\n";

	my $realign = "$gseq_l"."$gseq"."$gseq_r";
	my $bugle = "$tseq_l"."$rseq"."$tseq_r";
	$realign =~ tr/atgc/ATGC/;
	$bugle =~ tr/atgc/ATGC/;

	#print "$realign\n";
	#print "$bugle\n";
	#print "#########cccdddd#########\n";

	my $chr_id = $coor[0];
	my $r_start = $coor[1];
	my $r_end = $coor[2];
	my $gap_number = 0;
	my $mis = 0;

	($mis, $gap_number) = align($realign,$bugle);
	if($gap_number <= $gap_number_cutoff && $mis <= $gap_mismatch_cutoff){
		print "$chr_id\t$r_start\t$r_end\t$gap_number\t$mis\t$grna\t$realign\t$ali_str\t$bugle\t$rid\t-\n";
	}
}
close IN;

#################################
sub reverse_comp{
	my $seq = shift;	
	$seq =~ tr/atgc/ATGC/;
	$seq =~ tr/ATGC/TACG/;
	$seq = reverse($seq);
	return $seq;
}

################################
sub align{
        my $seq1 = shift;
        my $seq2 = shift;
        my @s1 = split //,$seq1;
        my @s2 = split //,$seq2;
        my $mismatch = 0;
	my $gap_count = 0;
        for(my $k = 0; $k < length($seq1); $k++){
                if($s1[$k] ne $s2[$k]){
			if($s1[$k] eq "N"){
				next;
			}elsif($s1[$k] eq "-" or $s2[$k] eq "-"){
				$gap_count++;
			}else{
                        	$mismatch++;
			}
                }
        }
        return $mismatch,$gap_count;
}

#################################
sub getSeq{
	my $c = shift;
	my $h = shift;
	
	my $chr_id = $c->[0];
	my $s = $c->[1];
	my $e = $c->[2];
	my $seq = $h->{$chr_id};

	my $target_seq = substr($seq,$s-1, $e-$s+1);
	return $target_seq;
}
