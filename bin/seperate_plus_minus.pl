#!/usr/bin/perl -w

open IN,"$ARGV[0]" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	if($t[5] eq "-"){
		print STDERR join "\t",@t;
		print STDERR "\n";
	}elsif($t[5] eq "+"){
		print join "\t",@t;
		print "\n";
	}
}
