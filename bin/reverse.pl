#!/usr/bin/perl -w

use strict;

my $index2 = "GCGATCTA";

open IN,"$ARGV[0]" || die $!;
while(<IN>){
	chomp;
	my $id = $_;
	my $seq = <IN>;
	my $id2 = <IN>;
	my $qual = <IN>;
	chomp($seq);
	chomp($id2);
	chomp($qual);
	#my $count = ($seq =~ m/$index2/g);
	#print "$count\n";
	#next if($count > 1);
	while($seq =~ m/$index2/g){
		my $len = length($index2);
		my $end = pos($seq);
		my $end2 = $end - 16 - 33 - 1;
		my $start = 27;

		my $subseq = substr($seq, $start, $end2-$start);
		my $subqual = substr($qual, $start, $end2-$start);
		print "$id\n$subseq\n$id2\n$subqual\n" if(length($subseq) > 0);
	}
}
close IN;
