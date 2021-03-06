#!/usr/bin/perl -w

use strict;

#id      read_count      no_bulge        bulge   target  target_realign  chrom   start   end
#chr9:110103647-110103748,chr9:110103677-110103778       1887    GGGAAAGACCCAGCATCCGTGGG         GGGAAAGACCCAGCATCCGTNGG         chr9    110103687       110103710

my %site;

my $i = 0;
foreach my $f_str (@ARGV){
	my ($f, $sample_id) = split /:/,$f_str;
	open IN,"$f" || die $!;
	<IN>;
	while(<IN>){
		chomp;
		my @t = split /\t/;
		my $id = "$t[6]:$t[7]-$t[8]";
		$site{$id}{count} = [];
		$site{$id}{id} = [];
		$site{$id}{fa_id} = [];
		#push @{$site{$id}{count}},$t[4];
		#push @{$site{$id}{id}},$sample_id;
		#push @{$site{$id}{fa_id}},$t[3];
		$site{$id}{ref_count} = $t[1] if($i == 0);

		$site{$id}{no_bulge} = $t[2];
		$site{$id}{bulge} = $t[3];
		$site{$id}{target} = $t[4];
		$site{$id}{target_realign} = $t[5];
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

#use Data::Dumper;
#print Dumper(\%site);

my $j = 0;
foreach my $f_str (@ARGV){
	my ($f, $sample_id) = split /:/,$f_str;
	open IN,"$f" || die $!;
	<IN>;
	while(<IN>){
		chomp;
		my @t = split /\t/;
		my $id = "$t[6]:$t[7]-$t[8]";
		my $c = $site{$id}{count};
		my $id2 = $site{$id}{id};
		my $fa_id = $site{$id}{fa_id};
		
		$c->[$j] = $t[1];
		$id2->[$j] = $sample_id;
		$fa_id->[$j] = $t[0];
	}
	close IN;
	$j++;
}

#use Data::Dumper;
#print Dumper(\%site);
#use List::Util qw[min max];

#id      read_count      no_bulge        bulge   target  target_realign  chrom   start   end
#chr9:110103647-110103748,chr9:110103677-110103778       1887    GGGAAAGACCCAGCATCCGTGGG         GGGAAAGACCCAGCATCCGTNGG         chr9    110103687       110103710
print "id\tread_count\tno_bulge\tbulge\ttarget\ttarget_realign\tchrom\tstart\tend\tsample_id\tsample_count\n";
foreach my $id (keys %site){
	#print "$id\n";
	my ($chr, $start, $end) = ($1,$2,$3) if($id =~ /(\w+):(\d+)-(\d+)/);
	my $count_str = join ";",@{$site{$id}{count}};
	my $id_str = join ";",@{$site{$id}{id}};
	my $fa_id_str = join ";",@{$site{$id}{fa_id}};

	my $no_bulge = $site{$id}{no_bulge};
	my $bulge = $site{$id}{bulge};
	my $target = $site{$id}{target};
	my $target_realign = $site{$id}{target_realign};

	$site{$id}{ref_count} = 0 unless(defined $site{$id}{ref_count});
	print "$fa_id_str\t$site{$id}{ref_count}\t$no_bulge\t$bulge\t$target\t$target_realign\t$chr\t$start\t$end\t$id_str\t$count_str\n";
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
