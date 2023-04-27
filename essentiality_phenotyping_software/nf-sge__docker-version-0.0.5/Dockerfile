FROM ubuntu:18.04 

LABEL maintainer="vo1@sanger.ac.uk" \
      version="0.0.5" \
      description="nf-sge container"

MAINTAINER  Victoria Offord <vo1@sanger.ac.uk>

USER root

# Locale
ENV LC_ALL C
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

# ALL tool versions used by opt-build.sh
ENV VER_BCL2FASTQ="2-20-0"
ENV VER_FASTQC="0.11.8"
ENV VER_SEQPREP="1.2"
ENV VER_EMBOSS="6.6.0"
ENV VER_SEQKIT="0.11.0"
ENV VER_MULTIQC="1.7"
ENV VER_CUTADAPT="2.5"

RUN apt-get -yq update

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip

RUN apt-get install -yq --no-install-recommends \
curl \
unzip \
make \
gcc \
g++ \
zlib1g-dev \
libfindbin-libs-perl \
openjdk-8-jre \
libx11-dev 

RUN apt-get -yq update

RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm Data::Dumper
RUN cpanm Getopt::Long
RUN cpanm File::Basename
RUN cpanm -fn XML::DOM::XPath
RUN cpanm -fn Bio::SeqIO

ENV OPT /opt/wsi-t113
ENV PATH $OPT/bin:$OPT/FastQC:$PATH
ENV LD_LIBRARY_PATH $OPT/lib
ENV PERL5LIB $OPT/lib/perl5

ENV DISPLAY=:0

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

RUN mkdir -p $OPT/bin

ADD build/opt-build.sh build/
RUN bash build/opt-build.sh $OPT

RUN chmod a+rx $OPT/bin

USER ubuntu

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

WORKDIR /home/ubuntu
