#!/bin/bash

CONSERTING_ROOT_DIR=/xchip/gistic/Jeremiah/Projects/Jabba/CONSERTING/code

if [[ ! -e $CONSERTING_ROOT_DIR ]]; then
    echo "Change the CONSERTING_ROOT_DIR to the path on your system"
    exit 1
fi

## modify paths
export PATH=$PATH:$CONSERTING_ROOT_DIR
export CLASSPATH=$CLASSPATH:$CONSERTING_ROOT_DIR:$CONSERTING_ROOT_DIR/bambino_bundle_1.06.jar:$CONSERTING_ROOT_DIR/sam-1.65.jar:$CONSERTING_ROOT_DIR/mysql-connector-java-5.1.39.tar.gz

TUM=$1
NORM=$2

## check that files exist, and link them (stupid bai vs bam.bai issues
linkr $TUM tum.bam
linkr $NORM norm.bam

## convert BAM to WIG
time sh bam2wig.sh tum.bam tum.wig
time sh bam2wig.sh norm.bam norm.wig

#### STILL NEED TO FIGURE OUT WHAT IS NEXT






## function for symlinking BAMs
function linkr() {
   
    if [[ ! -e $1 ]]; then
	echo "BAM file: $1 does not exists. Exiting"
	exit 1
    fi

    ## find BAM index, whether .bai or .bam.bai
    BAI=`echo ${1} | sed s/bam$/bai/`
    if [[ ! -e $BAI ]]; then
	BAI="${TUM}.bai"
    fi
    
    ## check that BAM index exists
    if [[ ! -e $BAI ]]; then
	echo "BAM index file: $BAI does not exist. Exiting"
	exit 1
    fi

    ## symlink with naming convention required for conserting
    ln -s ${1} ${2} 
    ln -s $BAI ${2}.bai
} 
