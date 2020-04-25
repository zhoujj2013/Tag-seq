#!/usr/bin/perl -w

use strict;

my $index_f = $ARGV[0];
my $reads_f = $ARGV[1];
my $water_result4visual_f = $ARGV[2];

# fasta id to id index
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

# get count information from bed file
my %count;
open IN,"$reads_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	
	#my $id = "$t[0]:$t[1]";
	my $id = "$t[0]:$t[1]-$t[2]";
	$count{$id} = [$t[4],$t[5],$t[6],$t[7]];
}
close IN;


# the water_result4visual_f format
#grna    142     27      96.3%   96.3%   0.0     126.0   TTTAGGATTCCCGAGTAGCAGATGACC     TTCAGGATTCCCGAGTAGCAGATGACC     chr3    46399672        46399698        +
#grna    144     27      100.0%  100.0%  0.0     135.0   TTTAGGATTCCCGAGTAGCAGATGACC     TTTAGGATTCCCGAGTAGCAGATGACC     chr3    46415023        46415049        +
#
print "chr\tstart\tend\tfa_id\tid\tplus_plus\tplus_minus\tminus_plus\tminus_minus\t\t\tbi.sum.mi\t\t\t\t\t\t\t\t\t\tOff-Target Sequence\t\t\t\t\t\t\t\t\t\t\tTarget Sequence\tstrand\n";
open IN,"$water_result4visual_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $off=$t[8];
	my $target = $t[7];
	my $id = $t[1];
	my $fa_id = $index{$id};
	my ($chr,$start,$end,$strand) = ($t[-4],$t[-3],$t[-2],$t[-1]); #if($fa_id =~ /(\w+):(\d+)-(\d+)/);
	my $c = $count{$index{$id}};
	my $sum = $c->[0] + $c->[1] + $c->[2] + $c->[3];
	print "$chr\t$start\t$end\t$fa_id\t$id\t$c->[0]\t$c->[1]\t$c->[2]\t$c->[3]\t\t\t$sum\t\t\t\t\t\t\t\t\t\t$off\t\t\t\t\t\t\t\t\t\t\t$target\t$strand\n";
}
close IN;

