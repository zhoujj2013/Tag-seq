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

        Remove ODN tag for read2 (allow mismatchs).
	But this script is very slow (1M, 10mins), I have obsoleted.
        Author: zhoujj2013\@gmail.com
        Usage: $0 xxx.r1.fq xxxx.r2.fq <ODN tag, ATCG> <max_mismatch> <outdir> <prefix>

USAGE
print "$usage";
exit(1);
};


my $r1_f = shift;
my $r2_f = shift;
my $odn = shift;
my $max_mismatch = shift;
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
	chomp($seq);
	chomp($qual);

	#print $id1;
	
	my $start = align($odn,$seq,$max_mismatch);
	#print "$start\n";
	if($start eq "NA"){
		$filtered{$id} = 1;
		next;
	}else{
		my $len = length($seq);
		my $new_seq = substr($seq,$start,$len-1);
		my $new_qual = substr($qual,$start,$len-1);
		if(length($new_seq) > 0){
			print OUT2 "$id1";
			print OUT2 "$new_seq\n";
			print OUT2 "$id2";
			print OUT2 "$new_qual\n";
		}else{
			$filtered{$id} = 1;
			next;
		}
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

###############
sub align {
	my $odn_seq = shift;
	my @o = split //,$odn_seq;

	my $q_seq = shift;
	my $max_mis = shift;

	#print "$odn_seq\t$q_seq\t$max_mis\n";
	
	my $odn_seq_len = length($odn_seq);
	my $q_seq_len = length($q_seq);

	my %k;
	for(my $i = 0; $i < $q_seq_len - $odn_seq_len + 1; $i++){
		my $s = $i+$odn_seq_len;
		$k{$s} =  substr($q_seq,$i,$odn_seq_len);
		#print length($k{$s})."\n";
	}
	#print scalar(keys %k)."\n";
	#print "#############\n";
	#print Dumper(\%k);

	my $pos = "NA";
	foreach my $key (keys %k){
		my @q = split //,$k{$key};
		my $mis = 0;
		for(my $j = 0; $j < scalar(@q); $j++){
			if($q[$j] ne $o[$j]){
				$mis++
			}
		}
		if($mis <= $max_mis){
			$pos = $key;
		}
	}
	return $pos;
}
