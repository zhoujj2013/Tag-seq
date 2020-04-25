#!/usr/bin/perl -w

use strict;

open IN,"$ARGV[0]" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $strand = $t[5];
	my $s;
	my $e;
	if($strand eq "-"){
		$s = $t[2]-1;
		$e = $s + 1;
	}elsif($strand eq "+"){
		$s = $t[1];
		$e = $t[1] + 1;
	}
	$t[1] = $s;
	$t[2] = $e;
	print join "\t",@t;
	print "\n";
}
close IN;
