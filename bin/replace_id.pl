#!/usr/bin/perl -w

use strict;

my %id_index;
open IN,"$ARGV[0]" || die $!;
$/ = ">"; <IN>; $/ = "\n";
my $i = 1;
while(<IN>){
	chomp;
	my $id = $1 if(/^(\S+)/);
	$/ = ">";
	my $seq = <IN>;
	chomp($seq);
	$/ = "\n";
	$id_index{$i} = $id;
	print ">$i\n";
	print "$seq";
	$i++;
}
close IN;

foreach my $k (keys %id_index){
	print STDERR "$k\t$id_index{$k}\n";
}
