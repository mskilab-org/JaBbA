#!/bin/bash

## get the tarball
wget http://ftp.stjude.org/pub/software/conserting/CONSERTING_code.tgz
tar xzfv CONSERTING_code.tgz
cd code

## get the bigwig to wig script
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig
chmod +x bigWigToWig

## get the mappability track for 100mers
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig

## convert bigwig to wiggle
./bigWigToWig wgEncodeCrgMapabilityAlign100mer.bigWig wgEncodeCrgMapabilityAlign100mer.wiggle

## make the mapability track thing
if hash java 2>/dev/null; then
    java -Xmx8g -cp ./ ExtractMappability Mapability wgEncodeCrgMapabilityAlign100mer.wiggle hg19.sizes 100
else
    echo "Need to have Java in PATH"
    exit 1
fi

## get the assembly files
mkdir -p HG19_ASSEMBLY
cd HG19_ASSEMBLY
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*
gunzip *gz
for f in *fa; do mv "$f" $(echo "$f" | sed 's/^chr//g'); done ## rename to get rid of chr
cd ..

## construct the GC content track
mkdir -p HG19_GC
java -cp ./ FA2GC $SIZE HG19_ASSEMBLY HG19_GC

## download bambino
wget https://cgwb.nci.nih.gov/goldenPath/bambino/bambino_bundle_1.06.jar

## download sam-1.65
wget http://ftp.stjude.org/pub/software/conserting/sam-1.65.jar

## download mysql-connector-java-5.1.10-bin.jar
wget http://dev.mysql.com/get/Downloads/Connector-J/mysql-connector-java-5.1.39.tar.gz

## download latest GNU parallel with nify script
(wget -O - pi.dk/3 || curl pi.dk/3/ || fetch -o - http://pi.dk/3) | bash
