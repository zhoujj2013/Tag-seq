# Tag-seq

## Description of Tag-seq data analysis

For Tag-seq data analysis, we retained fragments that contain an intact Tag at the beginning of read2 (second of pair). Then, reads were mapped to the reference genome (hg19) using STAR2 after quality filtering, then PCR duplications were removed using UMI-tools. To identify candidate DSBs, the start mapping positions were grouped if the distance among them is less than ten bps, resulting in editing hotspots induced by RGNs. Then, the peaks with sufficient reads were detected in RGNs hotspot. Furthermore, the peaks with reads mapping to both + and - strands, or the same strand but amplified with both forward and reverse tag-specific primers, are flagged as sites of potential DSBs. The flanking regions of potential DSBs match gRNA identified as on-target sites using a Smith-Waterman local-alignment algorithm. Identified off-targets sorted by Tag-seq read count are annotated in a final output table and visualize as a pdf file.

## System requirements

Tag-seq runs under the Linux (i.e., Centos, see also https://www.centos.org/ for further details) on a 64-bit machine with at least 32 GB RAM.

Tag-seq requires PERL v5, R, Python 2.7, [pip](https://bootstrap.pypa.io/get-pip.py) and several python packages listed in [python.package.requirement.txt](https://github.com/zhoujj2013/Tag-seq/blob/master/python.package.requirement.txt).;

Tag-seq also requires some third-party packages:

[STAR aligner](https://github.com/alexdobin/STAR)  
[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
[AdapterRemoval](https://github.com/MikkelSchubert/adapterremoval)  
[BEDTOOLS](https://bedtools.readthedocs.io/en/latest/)  
[SAMTOOLS](http://samtools.sourceforge.net/)  
[PICARD](https://broadinstitute.github.io/picard/)  
[umi_tools](https://github.com/CGATOxford/UMI-tools)  
[bedops](https://bedops.readthedocs.io/en/latest/)  
water in [EMBOSS](http://emboss.sourceforge.net/download/)  
[RIdeogram](https://github.com/TickingClock1992/RIdeogram)  

Tag-seq have been tested in CentOS release 7.4 (Linux OS 64 bit).

## Run by Docker

### Install docker

```
# update source
sudo apt-get update

# apt through https
sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
    software-properties-common

# add GPG key
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

# get the stable version
sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"

# update source
sudo apt-get update

# install docker-ce
sudo apt-get install docker-ce

# set user information, so that we can run docker without sudo
sudo usermod -a -G docker $USER

# exit and login again
```

### create docker image

```
git clone https://github.com/zhoujj2013/Tag-seq.git --depth 1
cd Tag-seq/docker/
docker build -t ubuntu:tagseq .
docker run -i -t ubuntu:tagseq echo "hello world!"
cd -
```

### run Tag-seq

Prepare reference
```
# prepare reference
mkdir ref && cd ref
rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz .
gunzip hg19.fa.gz
samtools faidx hg19.fa
rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes .
mkdir star_index
/path_to_star/STAR  --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ../hg19.fa --runThreadN 8
cd ../../
```

Prepare fastq dataset
```
# prepare fastq dataset as follow:
$tree data/
data/
├── 4P_L_R1.fq
├── 4P_L_R2.fq
├── 4P_R_R1.fq
├── 4P_R_R2.fq
└── README.txt

0 directories, 5 files
```

Run Tag-seq
```
cd Tag-seq/test
gunzip data/*.fq.gz
docker run -v /path_to/Tag-seq/test:/mnt/tagseq -v /path_to/ref:/mnt/tagseq/ref -w /mnt/tagseq -i -t ubuntu:tagseq perl /docker_main/software/Tag-seq/bin/run_guideseq.pl ./config.docker.txt all
#check the result
```

## Installation

### Get Tag-seq pipeline
```
git clone https://github.com/zhoujj2013/Tag-seq.git --depth 1
```


### Install require python packages

```
cd ./Tag-seq
pip install -r python.package.requirement.txt --user
```

### Preparation

Download reference genome and build index.

```
# download genome
mkdir hg19
cd hg19
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes
gunzip hg19.fa.gz

# build genome index
/path_to/STAR  --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ./hg19.fa --runThreadN 32

```

### Run test

If you have obtained the reference genome, STAR index, you can run test to examine whether the package works well (the test dataset is placed in ./test directory within Tag-seq).

Tag-seq requires a [sgrna.lst](https://github.com/zhoujj2013/Tag-seq/blob/master/test/sgrna.lst) and a configure file containing paths of input files, sgRNA, Tag primers and genome etc. (See [config.TEST.txt](https://github.com/zhoujj2013/Tag-seq/blob/master/test/config.TEST.txt) for more details.)

```
cd test

# run
########## the content of work.sh #########
# gunzip data/*.fq.gz
# perl ../bin/run_guideseq.pl config.TEST.txt all > config.TEST.log 2>config.TEST.err
###########################################
sh work.sh

# around 30 mins.
# you can check the report in out.XXX/sgrna_id.find.target/.
# you can identify off-targets for multiple sgRNA simultaneously.

```
### Result

#### 1. QC statistics

you can check [stat.txt](https://github.com/zhoujj2013/Tag-seq/blob/master/stat.txt).

#### 2. Information of potential targets in bed format

```
chr1    10111   10112   AAVS1.E_minus_minus_2_9,AAVS1.E_plus_minus_1_13 0       29      0       12
chr1    55903742        55903743        AAVS1.E_minus_minus_3669_8,AAVS1.E_plus_minus_5324_6    0       11      0       17
chr1    68164302        68164303        AAVS1.E_minus_minus_4802_6,AAVS1.E_plus_plus_6944_6     8       0       0       7
chr1    111700139       111700140       AAVS1.E_minus_minus_7763_6,AAVS1.E_plus_plus_11377_6    9       0       0       5
chr1    121478642       121478643       AAVS1.E_minus_plus_8420_6,AAVS1.E_plus_plus_12435_6     9       0       7       0
```

Column 1: chromosome  
Column 2: start  
Column 3: end  
Column 4: id  
Column 5: read count for plus strand in plus library  
Column 6: read count for minus strand in plus library  
Column 7: read count for plus strand in minus library  
Column 8: read count for minus strand in minus library  

#### 3. Potential off-targets

Illustrate of off-targets sites and read count.

![off-targets](https://upload-images.jianshu.io/upload_images/4180410-e4d77af6060e2f8a.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

#### 4. Read counts across sgRNA in target and off-target sites

![sites](https://upload-images.jianshu.io/upload_images/4180410-35ecdde6cda3c9c4.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240) 

#### 5. Global view of target and off-target sites

![global](https://upload-images.jianshu.io/upload_images/4180410-5303467568095e98.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


## Tag-seq Runtime

The running time of Tag-seq depends on the size of sequencing depth (For 30M flagments, it takes 30mins). 

## Please cite

1. xxxx Tag-seq (submitted)
