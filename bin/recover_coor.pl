#!/usr/bin/perl -w

use strict;

my $water_f = shift;
my $bed_f = shift;
my $index_f = shift;

my %index;
open IN,"$index_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	$index{$t[0]} = $t[1];
}
close IN;


my %bed;
open IN,"$bed_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = "$t[0]:$t[1]-$t[2]";
	$bed{$id} = \@t;
}
close IN;


open IN,"$water_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $grna = shift @t;
	my $tar = shift @t;

	my $bed_line = $bed{$index{$tar}};
	push @{$bed_line},@t;
	print join "\t",@{$bed_line};
	print "\n";
}
close IN;


