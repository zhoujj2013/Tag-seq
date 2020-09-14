#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);
use List::Util qw(min);

&usage if @ARGV<1;

sub usage {
        my $usage = << "USAGE";

        This script design for demultiplex of fastq files.
        Author: zhoujj2013\@gmail.com 
        Date: 2020.09.14
        Usage: $0 <prefix> <r1.fq> <r2.fq>

USAGE
print "$usage";
exit(1);
};

my @read_index;

my $prefix = shift;
my $r1_f = shift;
my $r2_f = shift;

my %barcode;
foreach my $f (@ARGV){
	my $fn1 = "$prefix.$f"."_r1";
	open "$fn1",">","$prefix.$f.r1.fq" || die $!;
	my $fn2 = "$prefix.$f"."_r2";
	open "$fn2",">","$prefix.$f.r2.fq" || die $!;
	$barcode{$f} = ();
}

open UN1,">","$prefix.undetermine.r1.fq" || die $!;
open UN2,">","$prefix.undetermine.r2.fq" || die $!;

#use Data::Dumper;
#print STDERR Dumper(\%barcode);

my %i2barcode;
my $i = 1;
open IN,"$r1_f" || die $!;
while(<IN>){
	chomp;
	my $id1 = $_;
	my $seq = <IN>;
	my $id2 = <IN>;
	my $qual = <IN>;
	
	my $index = substr($seq,8,8);
	
	my $fn = "$prefix.$index"."_r1";

	if(exists $barcode{$index}){
		#print $fn,"$id1\n";
		#push @{$barcode{$index}}$i;
		$i2barcode{$i} = $index;
		$fn->print("$id1\n"."$seq"."$id2"."$qual");
		#print $fn, "$id1\n"."$seq"."$id2"."$qual";
	}else{
		my $min = 8;
		my $sample_id = "";
		foreach my $b (keys %barcode){
			#print STDERR "$b\t$index\n";
			my $mismat = align($index,$b);
			if($mismat < $min){
				$sample_id = $b;
			}
		}
		if($min > 1){
			print UN1 "$id1\n"."$seq"."$id2"."$qual";
		}else{
			$i2barcode{$i} = $index;
			$fn->print("$id1\n"."$seq"."$id2"."$qual");
		}
	
	}
	$i++;
}
close IN;

#use Data::Dumper;
#print STDERR Dumper(\%i2barcode);

my $j = 1;
open IN,"$r2_f" || die $!;
while(<IN>){
	chomp;
	my $id1 = $_;
	my $seq = <IN>;
	my $id2 = <IN>;
	my $qual = <IN>;
	if(exists $i2barcode{$j}){
		my $fn = "$prefix.$i2barcode{$j}"."_r2";
		$fn->print("$id1\n"."$seq"."$id2"."$qual");
	}else{
		print UN2 "$id1\n"."$seq"."$id2"."$qual";
	}
	$j++;
}
close IN;


foreach my $f (@ARGV){
	my $fn1 = "$f"."_r1";
	my $fn2 = "$f"."_r2";
	close $fn1;
	close $fn2;
}
close UN1;
close UN2;


sub align{
	my $seq1 = shift;
	my $seq2 = shift;
	my @s1 = split //,$seq1;
	my @s2 = split //,$seq2;
	my $mismatch = 0;
	for(my $k = 0; $k < length($seq1); $k++){
		if($s1[$k] ne $s2[$k]){
			$mismatch++;
		}
	}
	return $mismatch;
}

