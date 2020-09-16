#!/usr/bin/perl -w

use strict;

my $gap_raw_f = shift;
my $gapBed_f = shift;
my $nongapBed_f = shift;

my %raw;
open IN,"$gap_raw_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = "$t[0]:$t[1]-$t[2]";
	$raw{$id} = \@t;
}
close IN;

my %gapBed;
open IN,"$gapBed_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = $t[3];
	$gapBed{$id} = \@t;
}
close IN;

print "id\tread_count\tno_bulge\tbulge\ttarget\ttarget_realign\tchrom\tstart\tend\n";
my @output;

# for gap overlap with nongap entry
my %nongapBed;
open IN,"$nongapBed_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = $t[3];
	my @id = split /,/,$id;
	#chr9    110103687       110103710       chr9:110103647-110103748,chr9:110103677-110103778       1887    GGGAAAGACCCAGCATCCGTGGG GGGAAAGACCCAGCATCCGTGGG +
	my @out;
	foreach my $pos (@id){
		$nongapBed{$pos} = 1;
		if(exists $gapBed{$pos}){
			my $raw_line = $raw{$pos};
			#chr12   5555155 5555256 1       3       GGGAAAGACCCAGCATCCGTGGG GGGAAAGACCCAGCATCCGTGGG |.|||||| |||||||||.|.|| GAGAAAGA-CCAGCATCCATAGG 105     +
			my $bulge = $raw_line->[8];
			my $re_align = $raw_line->[6];
			# if overlap, we only retain the non-overlap one.
			#@out = ($id,$t[4],$t[5],$bulge,$t[6],$re_align);
			# id	read_count	no_bulge	bulge	target	target_realign	chrom	start	end
			@out = ($id,$t[4],$t[5],"",$t[6],"",$t[0],$t[1],$t[2]);
			last;
		}else{
			@out = ($id,$t[4],$t[5],"",$t[6],"",$t[0],$t[1],$t[2]);
			last;
		}
	}
	push @output,\@out;
}
close IN;

#use Data::Dumper;
#print Dumper(\@output);


# for gap not overlapped with nongap entry
foreach my $id (keys %raw){
	unless(exists $nongapBed{$id}){
		#chr8    66842977        66843078        1       5       GGGAAAGACCCAGCATCCGTGGG GGGAAAGACCCAGCATCCGTGGG ||.|.||||.|||| .||.|||| GGCACAGACTCAGC-CCCCTGGG 410       +
		my $l = $raw{$id};

		next unless(exists $gapBed{$id});

		#if(exists $gapBed{$id}){
		#	print STDERR "Yes!\n";
		#}else{
		#	print STDERR "No!\n";
		#}

		my $gl = $gapBed{$id};
		my $read_count = $gl->[4];
		my $no_bulge = "";
		my $bulge = $l->[8];
		my $target = $l->[5];
		my $re_align = $l->[6];
		#chr8:66842977-66843078                             410         GGCACAGACTCAGC-CCCCTGGG  GGGAAAGACCCAGCATCCGTGGG  GGGAAAGACCCAGCATCCGTGGG
		my @out = ($id,$read_count,$no_bulge,$bulge,$target,$re_align,$l->[0], $l->[1], $l->[2]);
		push @output,\@out;
	}
}

my @output_sorted = sort {$b->[1] <=> $a->[1]} @output;

 
foreach my $l (@output_sorted){
	print join "\t",@$l;
	print "\n";
}
