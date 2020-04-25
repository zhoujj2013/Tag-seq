#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);

&usage if @ARGV<1;

sub usage {
        my $usage = << "USAGE";

        Prepare file for Rideogram.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 xx.merged.confirmed xx.visualization.offtarget.bed > xxx.sites.rideogram.txt

USAGE
print "$usage";
exit(1);
};

my $confirmed_f = shift;
my $offtarget_f = shift;

`bedtools intersect -a $confirmed_f -b $offtarget_f -u > ontarget.bed`;
`bedtools intersect -a $confirmed_f -b $offtarget_f -v > offtarget.bed`;

#Type    Shape   Chr     Start   End     color
#offtarget       triangle        chr1    10110   10111   ff7f00
#offtarget       triangle        chr1    121742261       121742262       ff7f00

print "Type\tShape\tChr\tStart\tEnd\tcolor\n";
open IN,"./ontarget.bed" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	print "Mapped_sgRNA\tbox\t$t[0]\t$t[1]\t$t[2]\t33a02c\n"
}
close IN;

open IN,"./offtarget.bed" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	print "Filtered\ttriangle\t$t[0]\t$t[1]\t$t[2]\tff7f00\n"
}
close IN;
