#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);

&usage if @ARGV<2;

sub usage {
        my $usage = << "USAGE";

        This script is designed for running Tag-seq analysis.
        Author: zhoujj2013\@gmail.com 2020/3/25
        Last update: 2020/8/31
        Usage: $0 <config.txt> <all|create_makefile|align|find_target>

USAGE
print "$usage";
exit(1);
};

my $conf=shift;
my $step=shift;

my %conf;
&load_conf($conf, \%conf);

$conf{OUTDIR} = abs_path($conf{OUTDIR});
$conf{FORWARD_LIB_R1} = abs_path($conf{FORWARD_LIB_R1});
$conf{FORWARD_LIB_R2} = abs_path($conf{FORWARD_LIB_R2});
$conf{REVERSE_LIB_R1} = abs_path($conf{REVERSE_LIB_R1});
$conf{REVERSE_LIB_R2} = abs_path($conf{REVERSE_LIB_R2});
$conf{CTRL} = abs_path($conf{CTRL});
$conf{BIN} = abs_path($Bin);
$conf{GRNA} = abs_path($conf{GRNA});

#print Dumper(\%conf);

mkdir $conf{OUTDIR} unless(-d "$conf{OUTDIR}");

if($step eq "all" || $step eq "create_makefile"){
# create sample.lst
open OUT,">","$conf{OUTDIR}/sample.lst" || die $!;
print OUT "$conf{PREFIX}_plus\t$conf{FORWARD_LIB_R1}\t$conf{FORWARD_LIB_R2}\tforward\t$conf{FORWARD_LIB_TAG}\n";
print OUT "$conf{PREFIX}_minus\t$conf{REVERSE_LIB_R1}\t$conf{REVERSE_LIB_R2}\treverse\t$conf{REVERSE_LIB_TAG}\n";
close OUT;


# create config file
open OUT,">","$conf{OUTDIR}/config.txt" || die $!;
print OUT "#sample\n";
print OUT "SAMPLE\t$conf{OUTDIR}/sample.lst\n";
print OUT "OUTDIR\t$conf{OUTDIR}\n";

print OUT "# parameter\n";
print OUT "MINLEN\t$conf{MINLEN}\n";
print OUT "READLEN\t$conf{READLEN}\n";
print OUT "MAXINS\t$conf{MAXINS}\n";
print OUT "THREAD\t$conf{THREAD}\n";
print OUT "ADAPTER\t$conf{ADAPTER}\n";

print OUT "#program\n";
print OUT "BIN\t$conf{BIN}\n";
print OUT "FASTQC\t$conf{FASTQC}\n";
print OUT "STAR\t$conf{STAR}\n";
print OUT "BEDTOOLS\t$conf{BEDTOOLS}\n";
print OUT "SAMTOOLS\t$conf{SAMTOOLS}\n";
print OUT "PICARD\t$conf{PICARD}\n";
print OUT "AdapterRemoval\t$conf{AdapterRemoval}\n";
print OUT "umi_tools\t$conf{umi_tools}\n";
print OUT "water\t$conf{water}\n";
print OUT "bedops\t$conf{bedops}\n";

print OUT "# species\n";
print OUT "INDEX\t$conf{INDEX}\n";
print OUT "GENOME\t$conf{GENOME}\n";
print OUT "REF\t$conf{REF}\n";
print OUT "CHROMSIZE\t$conf{CHROMSIZE}\n";
close OUT;

# run guide-seq
chdir $conf{OUTDIR};
`perl $conf{BIN}/guideseq.v5.pl config.txt > $conf{PREFIX}.log 2>$conf{PREFIX}.err`;
} # create makefile END. 

if($step eq "all" || $step eq "align"){
	chdir $conf{OUTDIR};
	`cd $conf{OUTDIR}/$conf{PREFIX}\_plus/ && make > log 2>err && cd -`;
	`cd $conf{OUTDIR}/$conf{PREFIX}\_minus/ && make > log 2>err && cd -`;
	
	# get stat
	`perl $conf{BIN}/getMatrics.pl $conf{PREFIX}\_plus $conf{PREFIX}\_minus > stat.txt`;
} # alignment END


if($step eq "all" || $step eq "find_target"){
	# detect off-targets for each sgRNAs
	open IN,"$conf{GRNA}" || die $!;
	while(<IN>){
		chomp;
		my @t = split /\t/;
		my $id = $t[0];
		my $grna_seq = $t[1];
		my $pam = $t[2];
		$conf{PAM} = $pam;
		
		mkdir "$conf{OUTDIR}/$id.find.target" unless(-d "$conf{OUTDIR}/$id.find.target");
		chdir "$conf{OUTDIR}/$id.find.target";

		if($conf{CTRL} eq "none"){
			`touch blank.bed`;
			$conf{CTRL} = "./blank.bed";
		}

		# create gRNA.fa
		open OUT,">","$conf{OUTDIR}/$id.find.target/grna.fa" || die $!;
		print OUT ">grna\n";
		print OUT "$grna_seq\n";
		close OUT;

		open OUT, ">", "run.sh" || die $!;
		print OUT "perl $conf{BIN}/find_targetsite.v2.pl ../$conf{PREFIX}\_plus/02potentialTargets/$conf{PREFIX}\_plus.plus.proximal.sorted ../$conf{PREFIX}\_plus/02potentialTargets/$conf{PREFIX}\_plus.minus.proximal.sorted ../$conf{PREFIX}\_minus/02potentialTargets/$conf{PREFIX}\_minus.plus.proximal.sorted ../$conf{PREFIX}\_minus/02potentialTargets/$conf{PREFIX}\_minus.minus.proximal.sorted $conf{CTRL} $conf{OUTDIR}/$id.find.target/grna.fa $conf{PAM} $conf{BIN}/../data/$conf{GENOME}.blacklist.bed $conf{REF} $conf{CHROMSIZE} $conf{PREFIX}\n";
		close OUT;
		
		`perl $conf{BIN}/find_targetsite.v2.pl ../$conf{PREFIX}\_plus/02potentialTargets/$conf{PREFIX}\_plus.plus.proximal.sorted ../$conf{PREFIX}\_plus/02potentialTargets/$conf{PREFIX}\_plus.minus.proximal.sorted ../$conf{PREFIX}\_minus/02potentialTargets/$conf{PREFIX}\_minus.plus.proximal.sorted ../$conf{PREFIX}\_minus/02potentialTargets/$conf{PREFIX}\_minus.minus.proximal.sorted $conf{CTRL} $conf{OUTDIR}/$id.find.target/grna.fa $conf{PAM} $conf{BIN}/../data/$conf{GENOME}.blacklist.bed $conf{REF} $conf{CHROMSIZE} $conf{PREFIX}`;
		
		# check sites
		mkdir "sites" unless(-d "sites");
		chdir "sites";
		`perl $conf{BIN}/monitoring_offtargets.pl ../$conf{PREFIX}.parsing_water_for_visualization.offtarget.bed ../../$conf{PREFIX}\_plus/03visual/$conf{PREFIX}\_plus.plus.bdg ../../$conf{PREFIX}\_plus/03visual/$conf{PREFIX}\_plus.minus.bdg ../../$conf{PREFIX}\_minus/03visual/$conf{PREFIX}\_minus.plus.bdg ../../$conf{PREFIX}\_minus/03visual/$conf{PREFIX}\_minus.minus.bdg > monitoring_offtargets.log 2>monitoring_offtargets.err`;
		chdir "..";
		
		# check global site distribution
		mkdir "global" unless(-d "global");
		chdir "global";
		if($conf{GENOME} eq "hg19"){
			`cp $conf{BIN}/RIdeogram/data/human_karyotype.hg19.txt ./karyotype.txt`;
			`cp $conf{BIN}/RIdeogram/data/human_genedensity.hg19.txt ./genedensity.txt`;
		}elsif($conf{GENOME} eq "hg38"){
			`cp $conf{BIN}/RIdeogram/data/human_karyotype.hg38.txt ./karyotype.txt`;
			`cp $conf{BIN}/RIdeogram/data/human_genedensity.hg38.txt ./genedensity.txt`;
		}
		`perl $conf{BIN}/RIdeogram/prepare_target_offtarget_for_Rideogram.pl ../$conf{PREFIX}.all.sites.merged.confirmed ../$conf{PREFIX}.parsing_water_for_visualization.offtarget.bed > ./offtargets.txt`;
		
		open UL,">","./run.sh" || die $!;
		print UL "ulimit -s 16384\n";
		print UL "Rscript $conf{BIN}/RIdeogram//run_RIdeogram.R offtargets.txt $conf{PREFIX} > RIdeogram.log 2>RIdeogram.err\n";
		close UL;
		`sh ./run.sh`;
		chdir "..";
	}
}# find target END

chdir "..";
chdir "..";

#########################
sub load_conf
{
    my $conf_file=shift;
    my $conf_hash=shift; #hash ref
    open CONF, $conf_file || die "$!";
    while(<CONF>)
    {
        chomp;
        next unless $_ =~ /\S+/;
        next if $_ =~ /^#/;
        warn "$_\n";
        my @F = split"\t", $_;  #key->value
        $conf_hash->{$F[0]} = $F[1];
    }
    close CONF;
}
