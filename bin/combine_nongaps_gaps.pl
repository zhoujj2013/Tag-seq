#!/usr/bin/perl -w

use strict;

my $nongap_f = shift;
my $gap_f = shift;

my @offtargets;

my %nongap;
open IN,"$nongap_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my @id = split /,/,$t[3];
	foreach my $id (@id){
		$nongap{$id} = \@t;
	}
	push @offtargets,\@t;
}
close IN;

my %gap;
open IN,"$gap_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = $t[3];
	unless(exists $nongap{$id}){
		$t[3] = "$t[3]-gap";
		push @offtargets,\@t;
	}
}
close IN;

my @offtargets_sorted = sort {$b->[4] <=> $a->[4]} @offtargets;
foreach my $l (@offtargets_sorted){
	print join "\t",@$l;
	print "\n";
}
