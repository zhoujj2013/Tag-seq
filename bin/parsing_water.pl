#!/usr/bin/perl -w

use strict;


my %index;
open IN,"$ARGV[1]" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	$index{$t[0]} = $t[1];
}
close IN;

open IN,"$ARGV[0]" || die $!;
$/ = "#======================================="; <IN>; $/ = "\n";
while(<IN>){
	$/ = "#=======================================";
	my $con = <IN>;
	chomp($con);
	my $align = <IN>;
	chomp($align);
	$/ = "\n";

	# recover coordinates
	#print "$align\n";
	my @ali = split /\n/,$align;
	#print "$ali[4]\n";
	my ($s,$e) = ($2,$3) if($ali[4] =~ /(\d+)\s+(\d+)\s+\S+\s+(\d+)/);
	#print "$s\t$e\n";

	my $seq1 = $1 if($con =~ /# 1: (\S+)/);
	my $seq2 = $1 if($con =~ /# 2: (\S+)/);
	
	my ($chr,$start,$end) = ($1, $2, $3) if($index{$seq2} =~ /(.*):(\d+)-(\d+)/);
	my $new_start = $start + $s;
	my $new_end = $start + $e;

	my $length = $1 if($con =~ /# Length: (\S+)/);
	my $identity = $1 if($con =~ /# Identity:.*\((\S+)\)/);
	my $similarity = $1 if($con =~ /# Similarity:.*\((\S+)\)/);

	my $gap = 0;
	if($con =~ /# Gaps:.*\((\S+)\)/){
		$gap = $1;
	}

	my $score = $1 if($con =~ /# Score: (\S+)/);
	print "$seq1\t$seq2\t$length\t$identity\t$similarity\t$gap\t$score\t$chr\t$new_start\t$new_end\n";
}
close IN;


#
# Aligned_sequences: 2
# 1: target_seq_1
# 2: 1
# Matrix: EDNAFULL
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 7
# Identity:       6/7 (85.7%)
# Similarity:     6/7 (85.7%)
# Gaps:           0/7 ( 0.0%)
# Score: 26.0
#
#

