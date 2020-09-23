FROM ubuntu:latest
MAINTAINER zhoujj2013 "zhoujj2013@163.com"

# update the source
RUN sed -i s@/archive.ubuntu.com/@/mirrors.aliyun.com/@g /etc/apt/sources.list
RUN apt-get clean
RUN apt-get update

# install python
RUN apt-get install python2.7 -y --force-yes
RUN ln -s /usr/bin/python2.7 /usr/bin/python
RUN apt-get install curl -y --force-yes && curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
RUN python get-pip.py

# install star
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get install -y --no-install-recommends --force-yes \
    build-essential \
    python-dev \
    bzip2 \
    cmake \
    default-jdk \
    git \
    libnss-sss \
    libtbb2 \
    libtbb-dev \
    ncurses-dev \
    tzdata \
    unzip \
    wget \
    zlib1g \
    zlib1g-dev

# install FASTQC
RUN apt-get install -y --force-yes --no-install-recommends \
    adapterremoval \
    bedtools \
    samtools \
    picard \
    bedops \
    emboss \
    librsvg2-bin
    
#    fastqc \

#ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
#dpkg-reconfigure -f noninteractive tzdata

# Star aligner 
ENV star_version 2.6.1d
WORKDIR /docker_main
#ADD https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz /usr/bin/
ADD http://120.78.148.57/tagseq/star/${star_version}.tar.gz /usr/bin/
RUN tar -xzf /usr/bin/${star_version}.tar.gz -C /usr/bin/
RUN cp /usr/bin/STAR-${star_version}/bin/Linux_x86_64/* /usr/local/bin

RUN cd /docker_main / && \
    apt-get autoremove -y && \
    apt-get autoclean -y  && \
    apt-get clean

#rm -rf ${star_version}.tar.gz && \
# Set default working path
WORKDIR /docker_main

# fastqc
RUN mkdir software && cd ./software \
    && wget http://120.78.148.57/tagseq/fastqc/fastqc_v0.11.9.zip \
    && unzip fastqc_v0.11.9.zip \
    && chmod 755 FastQC/fastqc
WORKDIR /docker_main

# picard.jar
RUN cd ./software \
    && wget http://120.78.148.57/tagseq/picard/picard.jar \
    && chmod 755 ./picard.jar

# download Tag-seq
RUN cd ./software \
    && git clone https://gitee.com/zhoujj2013/Tag-seq.git --depth 1 \
    && cd ./Tag-seq \
    && pip install -i https://pypi.mirrors.ustc.edu.cn/simple/ -r ./python.package.requirement.txt \
    && cd /
