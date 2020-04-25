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

        Normailze to piont position.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 xx.proximal.sorted <plus/minus>

USAGE
print "$usage";
exit(1);
};

my $f=shift;
my $strand=shift;

open IN,"$f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	if($strand eq "plus"){
		$t[2] = $t[1]+1;
	}elsif($strand eq "minus"){
		$t[1] = $t[2]-1;
	}
	print join "\t",@t;
	print "\n";
}
close IN;
