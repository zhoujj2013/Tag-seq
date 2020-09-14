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

        This script designed for offtarget analysis of Cas9 based editting.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 <DeIE output dir1> <DeIE output dir3> <DeIE output dir3> ...

USAGE
print "$usage";
exit(1);
};

my $PAM = shift @ARGV;

my @s;
my @sample_names;
foreach my $d (@ARGV){
	my $abs_d = abs_path($d);
	my $prefix = basename($abs_d);
	$abs_d = dirname($abs_d);
	#my @f = glob("$d/*");
	#my $prefix;
	#foreach my $f (@f){
	#		if($f =~ /.*\/(.*)_minus/){
	#		$prefix = $1;
	#		}
	#}
	push @sample_names,$prefix;
	my $arr = "$abs_d/$prefix.parsing_water_for_visualization.offtarget.combined.forVis.bed:$prefix";
	push @s,$arr;
}

#print Dumper(\@s);

my $s_str = join " ",@s;
my $sample_names_str = join ".",@sample_names;
my $title_str = join ",",@sample_names;

`perl $Bin/integrate_target_sites4compare.pl $s_str > $sample_names_str.bed`;
#print "perl $Bin/integrate_target_sites4compare.pl $s_str > $sample_names_str.bed\n";

`python $Bin/visualization4compare.py --identified_file $sample_names_str.bed --outfile ./$sample_names_str --title $title_str --PAM $PAM`;
`rsvg-convert -a -f pdf -o $sample_names_str.pdf $sample_names_str.svg`;

#usage: visualization4compare.py [-h] --identified_file IDENTIFIED_FILE
#                                --outfile OUTFILE --title TITLE [--PAM PAM]



#print Dumper(\@s);
#foreach my $l (@s){
#	`perl $Bin/integrate_target_sites4compare.pl $l->[1]:$l->[0]  ../../out.CCR5.LB/find.target/CCR5LB.parsing_water_for_visualization.offtarget.bed:CCR5LB > for_draw.bed$l->[1]:$l->[0]`;
#}
