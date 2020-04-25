#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
#use lib "/home/zhoujj/my_lib/pm";
#use bioinfo;

&usage if @ARGV<1;

#open IN,"" ||die "Can't open the file:$\n";
#open OUT,"" ||die "Can't open the file:$\n";

sub usage {
        my $usage = << "USAGE";

        Description of this script.
        Author: zhoujj2013\@gmail.com
        Usage: $0 <bed> <prefix>

USAGE
print "$usage";
exit(1);
};

my ($f, $prefix) = @ARGV;


#my $f_name = basename($f,".bed");

#`mv $f $f.bk`;

my $i = 1;
open IN,"$f" || die $!;
while(<IN>){
        if(/^#|track/){
                print $_;
                next;
        }
        chomp;
        my @t = split /\t/;
        if(defined $t[3]){
                $t[3] = "$prefix\_$i";
        }else{
                push @t,"$prefix\_$i";
        }
        print join "\t",@t;
        print "\n";
        $i++;
}
close IN;

