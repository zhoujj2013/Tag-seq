#!/usr/bin/perl -w

use strict;

my $index_f = $ARGV[0];
my $reads_f = $ARGV[1];
my $water_result4visual_f = $ARGV[2];


#4845    chr2:132586497-132586548
#8012    chr6:124296980-124297031
#3181    chr15:67200406-67200457

my %index;
open IN,"$index_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	$index{$t[0]} = $t[1];
}
close IN;

my %count;
open IN,"$reads_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	
	#my $id = "$t[0]:$t[1]";
	my $id = "$t[0]:$t[1]-$t[2]";
	$count{$id} = $t[3];
}
close IN;

#grna    2484    27      66.7%   66.7%   0.0     54.0    CCTGCTGTTGGGTCTGATCACCCCAAA     CCTGCCCCTGAGTCTTAACCCACCATA
#grna    4736    27      96.3%   96.3%   0.0     126.0   CCTGCTGTTGGGTCTGATCACCCCAAA     CCTGGTGTTGGGTCTGATCACCCCAAA
#grna    6171    27      100.0%  100.0%  0.0     135.0   CCTGCTGTTGGGTCTGATCACCCCAAA     CCTGCTGTTGGGTCTGATCACCCCAAA


print "\t\t\t\t\t\t\t\t\t\t\tbi.sum.mi\t\t\t\t\t\t\t\t\t\tOff-Target Sequence\t\t\t\t\t\t\t\t\t\t\tTarget Sequence\t\n";
open IN,"$water_result4visual_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $off=$t[-1];
	my $target = $t[-2];
	my $id = $t[1];
	my $c = $count{$index{$id}};
	print "\t\t\t\t\t\t\t\t\t\t\t$c\t\t\t\t\t\t\t\t\t\t$off\t\t\t\t\t\t\t\t\t\t\t$target\t\n";
}
close IN;

