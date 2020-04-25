#!/usr/bin/perl -w

my $f = shift;
my $strand = shift;

open IN,"$f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my @n = split /_/,$t[3];
	my @new_name = (@n[0..1],$strand,@n[2..3]);
	$t[3] = join "_",@new_name;
	print join "\t",@t;
	print "\n";
}
close IN;

