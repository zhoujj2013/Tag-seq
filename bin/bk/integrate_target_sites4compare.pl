#!/usr/bin/perl -w

use strict;

#chr3    46415022        46415049        chr3:46414986-46415087  2389    TTTAGGATTCCCGAGTAGCAGATGACC     TTTAGGATTCCCGAGTAGCAGATGACC
#chr3    46399671        46399698        chr3:46399625-46399726  720     TTCAGGATTCCCGAGTAGCAGATGACC     TTTAGGATTCCCGAGTAGCAGATGACC

my %site;

my $i = 0;
foreach my $f_str (@ARGV){
	my ($f, $sample_id) = split /:/,$f_str;
	open IN,"$f" || die $!;
	while(<IN>){
		chomp;
		my @t = split /\t/;
		my $id = "$t[0]:$t[1]-$t[2]";
		$site{$id}{count} = [];
		$site{$id}{id} = [];
		$site{$id}{fa_id} = [];
		#push @{$site{$id}{count}},$t[4];
		#push @{$site{$id}{id}},$sample_id;
		#push @{$site{$id}{fa_id}},$t[3];
		$site{$id}{ref_count} = $t[4] if($i == 0);
		$site{$id}{off} = $t[5];
		$site{$id}{target} = $t[6];
	}
	close IN;
	$i++;
}

foreach my $f_str (@ARGV){
	foreach my $id (keys %site){
		push @{$site{$id}{count}},0;
		push @{$site{$id}{id}},"-";
		push @{$site{$id}{fa_id}},"-";
	}
}

my $j = 0;
foreach my $f_str (@ARGV){
	my ($f, $sample_id) = split /:/,$f_str;
	open IN,"$f" || die $!;
	while(<IN>){
		chomp;
		my @t = split /\t/;
		my $id = "$t[0]:$t[1]-$t[2]";
		my $c = $site{$id}{count};
		my $id2 = $site{$id}{id};
		my $fa_id = $site{$id}{fa_id};
		
		$c->[$j] = $t[4];
		$id2->[$j] = $sample_id;
		$fa_id->[$j] = $t[3];
	}
	close IN;
	$j++;
}

#use Data::Dumper;
#print Dumper(\%site);
#use List::Util qw[min max];
print "chr\tstart\tend\tfa_id\tcount\toff\ttarget\tid_str\tcount_str\n";
foreach my $id (keys %site){
	my ($chr, $start, $end) = ($1,$2,$3) if($id =~ /(\w+):(\d+)-(\d+)/);
	my $count_str = join ";",@{$site{$id}{count}};
	my $id_str = join ";",@{$site{$id}{id}};
	my $fa_id_str = join ";",@{$site{$id}{fa_id}};
	my $off = $site{$id}{off};
	my $target = $site{$id}{target};
	$site{$id}{ref_count} = 0 unless(defined $site{$id}{ref_count});
	print "$chr\t$start\t$end\t$fa_id_str\t$site{$id}{ref_count}\t$off\t$target\t$id_str\t$count_str\n";
}


#foreach my $f_str (@ARGV){
#	foreach my $id (keys %site){
#		push @{$site{$id}{count}},0;
#		push @{$site{$id}{id}},"-";
#	}
#}
#
##use Data::Dumper;
##print Dumper(\%site);
#
#foreach my $f_str (@ARGV){
#	my ($f, $sample_id) = split /:/,$f_str;
#	open IN,"$f" || die $!;
#	while(<IN>){
#		chomp;
#		my @t = split /\t/;
#		my $id = "$t[0]:$t[1]-$t[2]";
#		$site{$id}{count}
#		$site{$id}{id} = [];
#	}
#	close IN;
#}
#
