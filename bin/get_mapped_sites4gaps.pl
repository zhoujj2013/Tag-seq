#!/usr/bin/perl -w

use strict;

my $gap_bed_f = shift;
my $fw_fa = shift;
my $rv_fa = shift;
my $site_ext_bed = shift;

my %fw;
open IN,"$fw_fa" || die $!;
$/ = ">"; <IN>; $/ = "\n";
while(<IN>){
	my $id = $1 if(/^(\S+)/);
	$/ = ">";
	my $seq = <IN>;
	chomp($seq);
	$/ = "\n";
	$seq =~ s/\n//g;
	$fw{$id} = $seq;
}
close IN;


my %rv;
open IN,"$rv_fa" || die $!;
$/ = ">"; <IN>; $/ = "\n";
while(<IN>){
	my $id = $1 if(/^(\S+)/);
	$/ = ">";
	my $seq = <IN>;
	chomp($seq);
	$/ = "\n";
	$seq =~ s/\n//g;
	$rv{$id} = $seq;
}
close IN;

my %count;
open IN,"$site_ext_bed" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id_str = "$t[0]:$t[1]-$t[2]";
	my $c = $t[4]+$t[5]+$t[6]+$t[7];
	$count{$id_str} = $c;
}
close IN;

my %gap_record;
open IN,"$gap_bed_f" || die $!;
#chr2    219845025       219845126       1       3       GAGTCCGAGCAGAAGAAGAAGGG GAGTCCGAGCAGAAG-AAGAAGG |||.||||||||||| ||||.|| GAGGCCGAGCAGAAGAAAGACGG 89      +
#chr2    219845050       219845151       1       3       GAGTCCGAGCAGAAGAAGAAGGG GAGTCCGAGCAGAAG-AAGAAGG |||.||||||||||| ||||.|| GAGGCCGAGCAGAAGAAAGACGG 90      +
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $offtarget_seq = $t[8];
	my $offtarget_align = $t[8];
	$offtarget_seq =~ s/-//g;
	my $strand = $t[-1];
	my $id = "$t[0]:$t[1]-$t[2]";
	my $read_count = $count{$id};
	my $seq = "";
	if($strand eq "-"){
		$seq = $rv{$id};
	}elsif($strand eq "+"){
		$seq = $fw{$id};
	}
	#print "$id\n";
	#print "$offtarget_seq\n";
	#print "$seq\n";
	while($seq =~ /$offtarget_seq/g){
		#print "#######\n";
		my $len = length($offtarget_seq);
		my $end = pos($seq);
		my $start = $end - $len + 1;
		
		if($strand eq "+"){
			my $new_start = $t[1] + $start - 1;
			my $new_end = $new_start + $len;
			my @l = ($t[0],$new_start,$new_end,$id,$read_count,$t[5],$offtarget_align,$strand);
			my $new_id = "$t[0]:$new_start-$new_end";
			push @{$gap_record{$new_id}},\@l;
			#print "$t[0]\t$new_start\t$new_end\t$id\t$t[-2]\t$t[5]\t$offtarget_seq\t$strand\n";
		}elsif($strand eq "-"){
			my $new_start = $t[1] + (length($seq) - $end);
			my $new_end = $new_start + $len;
			my @l = ($t[0],$new_start,$new_end,$id,$read_count,$t[5],$offtarget_align,$strand);
			my $new_id = "$t[0]:$new_start-$new_end";
			push @{$gap_record{$new_id}},\@l;
			#print "$t[0]\t$new_start\t$new_end\t$id\t$t[-2]\t$t[5]\t$offtarget_seq\t$strand\n";
		}
	}
}
close IN;

#chr8    128801241       128801264       chr8:128801216-128801317        2861    GAGTCCTAGCAGGAGAAGAAGAG GAGTCCGAGCAGAAGAAGAAGGG +
#chr2    73160981        73161004        chr2:73160942-73161043  2784    GAGTCCGAGCAGAAGAAGAAGGG GAGTCCGAGCAGAAGAAGAAGGG +

#use Data::Dumper;
#print Dumper(\%gap_record);

# remove the duplicate sites
foreach my $p (sort keys %gap_record){
	my @l = sort {$b->[4] <=> $a->[4]} @{$gap_record{$p}};
	my $top_l = $l[0];
	print join "\t",@$top_l;
	print "\n";
	#foreach my $l (@l){
	#	print join "\t",@$l;
	#	print "\n";
	#}
}

