# Tag-seq

Data analysis pipeline for Tag-seq.

----------
## System requirements

Tag-seq runs under the Linux (i.e., Centos, see also https://www.centos.org/ for further details) on a 64-bit machine with at least 32 GB RAM.

Tag-seq requires Python 2.7, PERL v5, [pip](https://bootstrap.pypa.io/get-pip.py) and several python packages listed in [python.package.requirement.txt](https://github.com/zhoujj2013/Tag-seq/blob/master/python.package.requirement.txt).;

Tag-seq also requires STAR aligner, RIdeogram, FASTQC, AdapterRemoval, BEDTOOLS, SAMTOOLS, PICARD, umi_tools, water, bedops.
You must install them one by one.

Tag-seq have been tested in CentOS release 7.4 (Linux OS 64 bit).

## Installation

### Get Tag-seq pipeline
```
git clone https://github.com/zhoujj2013/Tag-seq.git --depth 1
```

### Preparation
Download reference genome and build index.

```
# download genome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz

# build genome index
/path_to/STAR  --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ./hg19.fa --runThreadN 32
```
### Run demo

If you have installed the Tag-seq package and obtained the supporting dataset, you can run demo to examine whether the package works well (the test dataset is placed in ./test directory within Tag-seq).

```
cd test
# create makefile
sh work.sh

# around 30 mins.
# you can check the report in out.XXX/find_targets/.

```

## Tag-seq Runtime

The running time of Tag-seq depends on the size of sequencing depth.

## Please cite

1. xxxx Tag-seq (underreview)


