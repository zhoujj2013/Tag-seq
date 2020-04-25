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

        Remove ODN tag for read2.
        Author: zhoujj2013\@gmail.com
        Usage: $0 xxxx.r2.fq <ODN tag, ATCG>

USAGE
print "$usage";
exit(1);
};


my $r1_f = shift;
my $r2_f = shift;
my $odn = shift;
my $outdir = shift;
my $prefix = shift;

open OUT1,">","$outdir/$prefix.rmODN.R1.fq" || die $!;
open OUT2,">","$outdir/$prefix.rmODN.R2.fq" || die $!;
open OUT3,">","$outdir/$prefix.rmODN.stat" || die $!;

my %filtered;

open IN,"$r2_f" || die $!;
while(<IN>){
	my $id1 = $_;
	my $id = $1 if(/(\S+)/); 
	my $seq = <IN>;
	my $id2 = <IN>;
	my $qual = <IN>;
	#print $id1;
	my $start = 0;
	while($seq =~ /$odn/g){
		$start = pos($seq);
	}
	chomp($seq);
	chomp($qual);
	my $len = length($seq);
	my $new_seq = substr($seq,$start,$len-1);
	my $new_qual = substr($qual,$start,$len-1);
	if($seq !~ /$odn/g){
		$filtered{$id} = 1;
		next;
	}elsif(length($new_seq) > 0){
		print OUT2 "$id1";
		print OUT2 "$new_seq\n";
		print OUT2 "$id2";
		print OUT2 "$new_qual\n";
	}else{
		$filtered{$id} = 1;
		next;
	}
}
close IN;

my $read_total = 0;
my $read_with_odn = 0;

open IN,"$r1_f" || die $!;
while(<IN>){
	$read_total = $read_total + 1;
	my $id1 = $_;
	my $id = $1 if(/(\S+)/); 
	my $seq = <IN>;
	my $id2 = <IN>;
	my $qual = <IN>;
	if(exists $filtered{$id}){
		next;
	}else{
		$read_with_odn = $read_with_odn + 1;
		print OUT1 "$id1";
		print OUT1 "$seq";
		print OUT1 "$id2";
		print OUT1 "$qual";
	}
}
close IN;

print OUT3 "Raw flagment count: $read_total\n";
print OUT3 "Read count with ODN: $read_with_odn\n";

close OUT1;
close OUT2;
close OUT3;

