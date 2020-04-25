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

        Design for get statistical matrics for Guide-seq analysis.
        Author: zhoujj2013\@gmail.com
        Usage: $0 sample_prefix1 sample_prefix2 sample_prefix3 ...

USAGE
print "$usage";
exit(1);
};


# for remove ODN
print "Retrive reads with ODN\n";
print "\tFiltered_Tag\tTotal_Tag\tPassed_Tag\tPassRate\n";
my @readsfilteredODN;
my @readstotalRaw;
my @readspassODN;
my @readsperODN;
foreach my $s (@ARGV){
	my @f = glob("$s/00datafilter/*.rmODN.stat");
	# format for *.rmODN.stat
	#
	#Raw flagment count: 5290357
	#Read count with ODN: 27218
	#
	my $f = $f[0];
	&checkExists($f);
	my $total = 0;
	my $pass = 0;
	open IN,"$f" || die $!;
	while(<IN>){
		chomp;	
		$total = $1 if(/Raw flagment count: (\d+)/);
		$pass = $1 if(/Read count with ODN: (\d+)/);
	}
	close IN;
	
	my $filtered = $total - $pass;
	my $per = sprintf("%.2f",($pass)*100/$total);
	
	push @readsfilteredODN,$filtered;
	push @readstotalRaw,$total;
	push @readspassODN,$pass;
	push @readsperODN,$per;	
}
ArrAsCol(\@ARGV,\@readsfilteredODN,\@readstotalRaw,\@readspassODN,\@readsperODN);
print "\n";

# for triming, all are fragment
print "Trim low quanlity and adaptor\n";
print "\tFiltered_Tag\tTotal_Tag\tPassed_Tag\tPassRate\n";
my @readsfiltered;
my @readstotal;
my @readspass;
my @readsper;
foreach my $s (@ARGV){
	my @f = glob("$s/00datafilter/*adapterRemoval.setting");
	my $f = $f[0];
	&checkExists($f);
	my $total = 0;
	open IN,"$f" || die $!;
	while(<IN>){
		chomp;
		next unless(/Total number of read pairs:/);
		$total = $1 if(/Total number of read pairs: (\S+)/);
	}
	close IN;
	
	my $pass = 0;
	my @f1 = glob("$s/01alignment/$s.Log.final.out");
	my $f1 = $f1[0];
	&checkExists($f1);
	open IN,"$f1" || die $!;
	while(<IN>){
		chomp;
		next unless(/Number of input reads/);
		#print $_,"\n";
		$pass = $1 if(/Number of input reads.*\s+(\d+)/);
	}
	close IN;
	
	#print $pass,"\n";

	my $filtered = $total - $pass;
	my $per = sprintf("%.2f",($pass)*100/$total);
	
	push @readsfiltered,$filtered;
	push @readstotal,$total;
	push @readspass,$pass;
	push @readsper,$per;	
}
ArrAsCol(\@ARGV,\@readsfiltered,\@readstotal,\@readspass,\@readsper);
print "\n";


#for mapping rate
print "Alignment Summary\n";
print "\tTotal_Tag\tAligned_Tag\tAlignmentRate\n";
my @aligned;
my @AlignRate;
my $i = 0;
foreach my $s (@ARGV){
	my $f = "$s/01alignment/$s.Log.final.out";
	&checkExists($f);
	my $readCleanTotal = $readspass[$i];
	my $rate;
	my $aligned;

	my $uniq_aligned=0;
	my $multiple_aligned=0;
	my $many_loci=0;
	my $total_sample=0;

	open IN,"$f" || die $!;
	while(<IN>){
		chomp;
		if(/Number of input reads.*\s+(\d+)/){
			$total_sample=$1;
		}
		if(/Uniquely mapped reads number.*\s+(\d+)/){
			$uniq_aligned=$1;
		}

		if(/Number of reads mapped to multiple loci.*\s+(\d+)/){
			$multiple_aligned=$1;
		}

		if(/Number of reads mapped to too many loci.*\s+(\d+)/){
			$many_loci=$1;
		}
		#$rate = $1 if(/(\S+)% overall alignment rate/);
		#$aligned = int($readCleanTotal * $rate * 0.01);
	}
	close IN;

	$aligned = $uniq_aligned + $multiple_aligned + $many_loci;
	$rate = sprintf("%.2f",($aligned/$total_sample)*100);
	push @AlignRate,$rate;
	push @aligned,$aligned;
	$i++;
}
$i=0;
ArrAsCol(\@ARGV,\@readspass,\@aligned,\@AlignRate);
print "\n";

#for deduce UMI
print "Duplication report (UMI)\n";
print "Duplication rate\n";
print "\tTotal_tag\tUniq_Tag\tDuplication rate(%)\n";
my @uniq_tag;
my @DupRate;
foreach my $s (@ARGV){
	my $f = "$s/01alignment/$s.Aligned.sortedByCoord.out.dedup.bam.flagstat";
	&checkExists($f);
	my $readMappedTotal = $aligned[$i];
	my $uniq;
	my $duprate;
	open IN,"$f" || die $!;
	while(<IN>){
		chomp;
		next unless(/paired in sequencing/);
		$uniq = $1 if(/^(\d+).*paired in sequencing/);
		$uniq = int($uniq/2);
		$duprate = sprintf "%.2f",(($readMappedTotal - $uniq)/$readMappedTotal)*100;
	}
	close IN;
	push @uniq_tag,$uniq;
	push @DupRate,$duprate;
	$i++;
}
$i=0;
ArrAsCol(\@ARGV,\@aligned,\@uniq_tag,\@DupRate);
print "\n";

#for filtering
print "Filtering report (mapping quality and tag filtering)\n";
print "Usable rate\n";
print "\tTotal_tag\tUsable_Tag\tFiltering rate(%)\n";
my @usable_tag;
my @filterRate;
foreach my $s (@ARGV){
	my $f = "$s/02potentialTargets/$s.Aligned.sortedByCoord.out.dedup.single.filtered.bam.flagstat";
	&checkExists($f);
	my $readMappedTotal = $uniq_tag[$i];
	my $usable;
	my $filterrate;
	open IN,"$f" || die $!;
	while(<IN>){
		chomp;
		next unless(/mapped \(/);
		$usable = $1 if(/^(\d+).*mapped \(/);
		$filterrate = sprintf "%.2f",(($readMappedTotal - $usable)/$readMappedTotal)*100;
	}
	close IN;
	push @usable_tag,$usable;
	push @filterRate,$filterrate;
	$i++;
}
$i=0;
ArrAsCol(\@ARGV,\@uniq_tag,\@usable_tag,\@filterRate);
print "\n";

print "Potential Target Num\n";
print "\tForward\tReverse\n";
my @filtered_cluster1;
my @filtered_cluster2;

foreach my $s (@ARGV){
	my $f1 = "$s/02potentialTargets/$s.plus.proximal.sorted";
	my $f2 = "$s/02potentialTargets/$s.minus.proximal.sorted";
	#my $raw_c = `wc -l $f | awk '{print \$1}'`;
	#chomp($raw_c);
	my $filtered_cluster1 = `awk '\$5 >= 5' $f1 | wc -l`;
        chomp($filtered_cluster1);
	push @filtered_cluster1,$filtered_cluster1;

	my $filtered_cluster2 = `awk '\$5 >= 5' $f2 | wc -l`;
        chomp($filtered_cluster2);
	push @filtered_cluster2,$filtered_cluster2;
}
ArrAsCol(\@ARGV,\@filtered_cluster1,\@filtered_cluster2);
print "\n";

sub ArrAsCol{
	my $arr = \@_;
	my $len = scalar(@{$arr->[0]});
	for(my $i = 0; $i < $len; $i++){
		my @l;
		foreach my $a (@$arr){
			push @l,$a->[$i];
		}
		print join "\t",@l;
		print "\n";
	}
}

sub checkExists{
	my $ffname = shift @_;
	unless(-e $ffname){
		print STDERR "$ffname not exists!\n";
		exit(1);
	}
}
