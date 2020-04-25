#!/usr/bin/perl -w

use strict;

# we optimise the primer design,
# so the ODN sequence is different from original one.
#
# After alignment, if the original sequence aligned to reverse strand, the sequence will be reversed and complemented in bam file. -- 2019.09
# so for plus library, we should seperate forward and reverse strand alignments.
# 1. 
 

my $tag = $ARGV[1]; # Input ODN: library sequence; for reverse alignment
my $type = $ARGV[2];

# set reverse ODN
# Reverse complement ODN, for forward alignment
my $rc_tag = $tag;
$rc_tag =~ tr/ATCG/TAGC/;
$rc_tag = reverse($rc_tag);
my $rc_tag_len = length($rc_tag);

print STDERR "$rc_tag\t$rc_tag_len\n";

#########################################################################
open IN,"$ARGV[0]" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $seq = $t[9];
	# reverse lib
	# GTTTAATTGAGTTGTCATATGTTAATAAC
	# 29bp
	# forward lib
	# ATACCGTTATTAACATATGACAACT
	# 25bp
	my $mis;
	my $mis_rc;
	

	# new scheme
	my $pos = $rc_tag_len*(-1);
	my $read_tag = substr($seq, $pos, $rc_tag_len); # reverse alignment
	my $read_tag_rc = substr($seq, 0, $rc_tag_len); # forward alignment
	
	if($type eq "forward"){
		$mis = align($read_tag, $tag); # reverse alignment
		$mis_rc = align($read_tag_rc, $rc_tag); # forward alignment
	}elsif($type eq "reverse"){
		$mis = align($read_tag, $tag); # reverse alignment
		$mis_rc = align($read_tag_rc, $rc_tag); # forward alignment
	}

	# old scheme2
	#if($type eq "forward"){
	#	#18bp
	#	my $rev_read_tag = substr($read_tag, 0, 18);
	#	my $rev_tag = substr($tag, 0, 18);
	#	$mis = align($rev_read_tag, $rev_tag);

	#	my $fw_read_tag_rc = substr($read_tag_rc, $pos, 18);
	#	my $fw_rc_tag = substr($rc_tag, $pos, 18);
	#	$mis_rc = align($fw_read_tag_rc, $fw_rc_tag);
	#}elsif($type eq "reverse"){
	#	#14bp
	#	my $rev_read_tag = substr($read_tag, 0, 14);
	#	my $rev_tag = substr($tag, 0, 14);
	#	$mis = align($rev_read_tag, $rev_tag);
	#	
	#	my $rev_read_tag_rc = substr($read_tag_rc, $pos, 14);
	#	my $rev_rc_tag = substr($rc_tag, $pos, 14);
	#	$mis_rc = align($rev_read_tag_rc, $rev_rc_tag);
	#}

	# old scheme
	#my $pos = $rc_tag_len*(-1);
	#my $read_tag = substr($seq, $pos, $rc_tag_len);
	#$mis = align($read_tag, $tag);
	#
	#my $read_tag_rc = substr($seq, 0, $rc_tag_len);
	#$mis_rc = align($read_tag_rc, $rc_tag);

	#if($type eq "forward"){
	#	my $read_tag = substr($seq, 0, 22);
	#	$mis = align($read_tag, $tag)
	#}elsif($type eq "reverse"){
	#	my $read_tag = substr($seq, 0, 22);
	#	$mis = align($read_tag, $tag)
	#}

	# setting missmatch
	if($mis > 0 && $mis_rc > 0){
		next;
	}else{
		print join "\t",@t;
		print "\n";
	}
}
close IN;


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
