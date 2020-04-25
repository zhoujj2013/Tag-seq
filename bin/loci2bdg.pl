#!/usr/bin/perl -w

use strict;

my %count;
open IN,"$ARGV[0]" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = "$t[0]#$t[1]#$t[2]";
	# seperate strand alignment
	my $strand = $t[-1];
	unless(exists $count{$id}){
		$count{$id}{"-"} = 0;
		$count{$id}{"+"} = 0;
	}
	if($strand eq "+"){
		$count{$id}{$strand}++;
	}else{
		$count{$id}{$strand}--;
	}
}
close IN;

foreach my $id (keys %count){
	my @id = split /#/,$id;
	my $plus = "+";
	my $minus = "-";
	print "$id[0]\t$id[1]\t$id[2]\t$count{$id}{$plus}\n";
	print STDERR "$id[0]\t$id[1]\t$id[2]\t$count{$id}{$minus}\n";
}
