#!/usr/bin/perl -w

use strict;

my $fw = shift; # waterman alignment result, *.water
my $rv = shift; # waterman alignment result, *.rev.water
my $fw_fa = shift; # input renamed fa, forward strand
my $rv_fa = shift; # input renamed fa, reverse strand
my $fw_renamed_index_f = shift;
my $rv_renamed_index_f = shift;
my $grna = shift; # grna fasta file
my $ref_length = shift; # grna length

my $pre_mismatch = 7;

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

# Read in sequence. why we need sequence??
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
# end

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
	
	my ($chr,$start,$end) = ($1, $2, $3) if($index_fw{$seq2} =~ /(.*):(\d+)-(\d+)/);
	my $new_start = $start + $s - 1;
	my $new_end = $start + $e;
	
	# check it, it is right.
	#print "$chr\t$new_start\t$new_end\t.\t.\t+\n";

	my $length = $1 if($con =~ /# Length: (\S+)/);
	my $identity = $1 if($con =~ /# Identity:.*\((\S+)\)/);
	my $similarity = $1 if($con =~ /# Similarity:.*\((\S+)\)/);

	# get gap information
	my $gap = 0;
	if($con =~ /# Gaps:.*\((.*)\)/){
		$gap = $1;
		$gap =~ s/%//g;
		$gap =~ s/ //g;
		$gap =~ s/\s+//g;
	}
	next if($gap > 0);
	next if($length/$ref_length < 0.5);
	if($ref_length > $length){
		# this part is very important
		# if the aligned sequence locate in marginal region, it is very difficult to process.
		#print "$seq1\t$seq2\n";
		my @sarr;
		my @seq = split //,$fwseq{$seq2};
		for(my $i = 0; $i < length($fwseq{$seq2}) - $ref_length + 1; $i++){
			my $subseq = substr($fwseq{$seq2},$i,$ref_length);
			my $mis = align($grna, $subseq);
			if($mis <= $pre_mismatch){
				#print "$seq1\t$seq2\t$mis\n";	
				my $score = $1 if($con =~ /# Score: (\S+)/);
				my $ref = $grna;
				my $offtarget = $subseq;
				$offtarget = uc($offtarget);
			
				# the coordinate aligned sgRNA
				$chr = $chr;
				my $new_new_start = $start + $i;
				my $new_new_end = $start + $i + $ref_length;	
				

				# this reverse complement will be elimilated
				#$ref = reverse($ref);
				#$ref =~ tr/TCGA/AGCT/;
				#$offtarget = reverse($offtarget);
				#$offtarget =~ tr/TCGA/AGCT/;

				print "$seq1\t$seq2\t$length\t$identity\t$similarity\t$gap\t$score\t$ref\t$offtarget\t$chr\t$new_new_start\t$new_new_end\t+\n";
			}
			#print STDERR substr($rvseq{$seq2},$i,$ref_length),"\n";
		}	
	}elsif($ref_length == $length){
		my $score = $1 if($con =~ /# Score: (\S+)/);
		my $ref = $1 if($con2 =~ /grna\s+\d+\s+(\w+)\s+\d+/);
		my $offtarget = $1 if($con2 =~ /\d+\s+\d+\s+(\w+)\s+\d+/);
		$offtarget = uc($offtarget);
		my $mis = align($grna, $offtarget);
		#print "$grna\t$offtarget\t$mis\n";
		if($mis <= $pre_mismatch){
			
			# this reverse complement will be elimilated
			#$ref = reverse($ref);
			#$ref =~ tr/TCGA/AGCT/;
			#$offtarget = reverse($offtarget);
			#$offtarget =~ tr/TCGA/AGCT/;
			print "$seq1\t$seq2\t$length\t$identity\t$similarity\t$gap\t$score\t$ref\t$offtarget\t$chr\t$new_start\t$new_end\t+\n";
		}
	}	
}
close IN;

# process water alignment (reverse part)
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
	if($con =~ /# Gaps:.*\((.*)\)/){
		$gap = $1;
		$gap =~ s/%//g;
		$gap =~ s/ //g;
		$gap =~ s/\s+//g;
	}
	next if($gap > 0);
	next if($length/$ref_length < 0.5);
	if($ref_length > $length){
		# this part is very important
		#print "$seq1\t$seq2\n";
		my @sarr;
		my @seq = split //,$rvseq{$seq2};
		for(my $i = 0; $i < length($rvseq{$seq2}) - $ref_length + 1; $i++){
			my $subseq = substr($rvseq{$seq2},$i,$ref_length);
			my $mis = align($grna, $subseq);
			if($mis <= $pre_mismatch){		
				#print "$seq1\t$seq2\t$mis\n";
				my $score = $1 if($con =~ /# Score: (\S+)/);
				my $ref = $grna;
				my $offtarget = $subseq;
				$offtarget = uc($offtarget);
				
				# the coordinate aligned sgRNA
				$chr = $chr;
				my $new_new_start = $end - $i - $ref_length;
				my $new_new_end = $end - $i;

				#
				#$ref = reverse($ref);
				#$ref =~ tr/TCGA/AGCT/;
				#$offtarget = reverse($offtarget);
				#$offtarget =~ tr/TCGA/AGCT/;
				#

				print "$seq1\t$seq2\t$length\t$identity\t$similarity\t$gap\t$score\t$ref\t$offtarget\t$chr\t$new_new_start\t$new_new_end\t-\n";
			}
			#print STDERR substr($rvseq{$seq2},$i,$ref_length),"\n";
		}
		
	}elsif($ref_length == $length){	
		my $score = $1 if($con =~ /# Score: (\S+)/);
		my $ref = $1 if($con2 =~ /grna\s+\d+\s+(\w+)\s+\d+/);
		my $offtarget = $1 if($con2 =~ /\d+\s+\d+\s+(\w+)\s+\d+/);
		$offtarget = uc($offtarget);
		my $mis = align($grna, $offtarget);
		if($mis <= $pre_mismatch){

			#$ref = reverse($ref);
			#$ref =~ tr/TCGA/AGCT/;
			#$offtarget = reverse($offtarget);
			#$offtarget =~ tr/TCGA/AGCT/;
			print "$seq1\t$seq2\t$length\t$identity\t$similarity\t$gap\t$score\t$ref\t$offtarget\t$chr\t$new_start\t$new_end\t-\n";
		}
	}	
}
close IN;
# end

#=======================================
#
# Aligned_sequences: 2
# 1: grna
# 2: 6167
# Matrix: EDNAFULL
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 26
# Identity:      25/26 (96.2%)
# Similarity:    25/26 (96.2%)
# Gaps:           0/26 ( 0.0%)
# Score: 121.0
#
#
#=======================================
#
#grna               1 CCTGCTGTTGGGTCTGATCACCCCAA     26
#                     ||||.|||||||||||||||||||||
#6167              42 cctggtgttgggtctgatcaccccaa     67
#
#
sub align{
        my $seq1 = shift;
        my $seq2 = shift;
        my @s1 = split //,$seq1;
        my @s2 = split //,$seq2;
        my $mismatch = 0;
        for(my $k = 0; $k < length($seq1); $k++){
                if($s1[$k] ne $s2[$k]){
                        $mismatch++;
                }
        }
        return $mismatch;
}
