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

        Create coverage figures for ontarget and offtarget sites.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 xx.visualization.offtarget.bed xx_plus.plus.bdg xx_plus.minus.bdg xx_minus.plus.bdg xx_minus.minus.bdg

USAGE
print "$usage";
exit(1);
};

my $offtarget_bed=shift;
my $plus_plus_bdg=shift;
my $plus_minus_bdg=shift;
my $minus_plus_bdg=shift;
my $minus_minus_bdg=shift;

`bedtools intersect -a $plus_plus_bdg -b $offtarget_bed -u > ./plus_plus.bdg`;
`bedtools intersect -a $plus_minus_bdg -b $offtarget_bed -u > ./plus_minus.bdg`;
`bedtools intersect -a $minus_plus_bdg -b $offtarget_bed -u > ./minus_plus.bdg`;
`bedtools intersect -a $minus_minus_bdg -b $offtarget_bed -u > ./minus_minus.bdg`;

my %pp;
open IN,"./plus_plus.bdg" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = "$t[0]#$t[1]#$t[2]";
	$pp{$id} = $t[3];
}
close IN;

my %pm;
open IN,"./plus_minus.bdg" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = "$t[0]#$t[1]#$t[2]";
	$pm{$id} = abs($t[3]);
}
close IN;

my %mp;
open IN,"./minus_plus.bdg" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = "$t[0]#$t[1]#$t[2]";
	$mp{$id} = $t[3];
}
close IN;

my %mm;
open IN,"./minus_minus.bdg" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $id = "$t[0]#$t[1]#$t[2]";
	$mm{$id} = abs($t[3]);
}
close IN;

open IN,"$offtarget_bed" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $chr = $t[0];
	my $s = $t[1];
	my $e = $t[2];
	# collection sequencing depth
	my $strand = $t[-1];
	my $offtarget_id = "$chr.$s.$e.$strand";
	open OUT,">","./$chr.$s.$e.$strand.wins.bed" || die $!;
	for(my $i = $s; $i < $e;$i++){
		my $interval_end = $i + 1;
		my $id = "$chr#$i#$interval_end";
		my $p_p = 0;
		my $p_m = 0;
		my $m_p = 0;
		my $m_m = 0;
		$p_p = $pp{$id} if(exists $pp{$id});
		$p_m = $pm{$id} if(exists $pm{$id});
		$m_p = $mp{$id} if(exists $mp{$id});
		$m_m = $mm{$id} if(exists $mm{$id});
		print OUT "$chr\t$i\t$interval_end\t$p_p\t$p_m\t$m_p\t$m_m\n";
	}
	close OUT;
	
	my $offtarget_seq = $t[5];
	my $ref_seq = $t[6];
	# use python to draw svg
	#print "python /home/zhoujj/project/02guide_seq/guideseq202003/visual_barchart.py $offtarget_seq $ref_seq $strand ./$chr.$s.$e.$strand.wins.bed $offtarget_id\n";
	`python $Bin/visual_barchart.py $offtarget_seq $ref_seq $strand ./$chr.$s.$e.$strand.wins.bed $offtarget_id`;
}
close IN;

