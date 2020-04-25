#!/usr/bin/perl -w

use strict;

my $raw_target_f = shift;

my %site;
open IN,"$raw_target_f" || die $!;
<IN>;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = "$t[0]:$t[1]-$t[2]";
	push @{$site{$id}},\@t;
}
close IN;

print "chr\tstart\tend\tfa_id\tid\tplus_plus\tplus_minus\tminus_plus\tminus_minus\t\t\tbi.sum.mi\t\t\t\t\t\t\t\t\t\tOff-Target Sequence\t\t\t\t\t\t\t\t\t\t\tTarget Sequence\tstrand\n";
foreach my $id (keys %site){
	my $lines = $site{$id};

	#print "$id\n";
	my ($chr,$start,$end) = ($1,$2,$3) if($id =~ /(\w+):(\d+)-(\d+)/);
	#print "$chr,$start,$end\n";
	my @fa_id;
	my @id_2;

	my $plus_plus = 0;
	my $plus_minus = 0;
	my $minus_plus = 0;
	my $minus_minus = 0;

	my $total_tags = 0;

	my $off;
	my $target; 
	my $strand;

	foreach my $l (@{$lines}){
		push @fa_id,$l->[3];
		push @id_2,$l->[4];
		$plus_plus = $plus_plus + $l->[5];
		$plus_minus = $plus_minus + $l->[6];
		$minus_plus = $minus_plus + $l->[7];
		$minus_minus = $minus_minus + $l->[8];
		$off = $l->[21];
		$target = $l->[32];
		$strand = $l->[33];
	}

	$total_tags = $plus_plus + $plus_minus + $minus_plus + $minus_minus;
	my $fa_id = join ",",@fa_id;
	my $id_2 = join ",",@id_2;
	print "$chr\t$start\t$end\t$fa_id\t$id_2\t$plus_plus\t$plus_minus\t$minus_plus\t$minus_minus\t\t\t$total_tags\t\t\t\t\t\t\t\t\t\t$off\t\t\t\t\t\t\t\t\t\t\t$target\t$strand\n";
}

#print "chr\tstart\tend\tfa_id\tid\tplus_plus\tplus_minus\tminus_plus\tminus_minus\t\t\tbi.sum.mi\t\t\t\t\t\t\t\t\t\tOff-Target Sequence\t\t\t\t\t\t\t\t\t\t\tTarget Sequence\t\n";
#open IN,"$water_result4visual_f" || die $!;
#while(<IN>){
#        chomp;
#        my @t = split /\t/;
#        my $off=$t[8];
#        my $target = $t[7];
#        my $id = $t[1];
#        my $fa_id = $index{$id};
#        my ($chr,$start,$end) = ($t[-4],$t[-3],$t[-2]); #if($fa_id =~ /(\w+):(\d+)-(\d+)/);
#        my $c = $count{$index{$id}};
#        my $sum = $c->[0] + $c->[1] + $c->[2] + $c->[3];
#        print "$chr\t$start\t$end\t$fa_id\t$id\t$c->[0]\t$c->[1]\t$c->[2]\t$c->[3]\t\t\t$sum\t\t\t\t\t\t\t\t\t\t$off\t\t\t\t\t\t\t\t\t\t\t$target\t\n";
#}
#close IN;
#
