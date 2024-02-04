#!/usr/bin/env sh

prct='0.05'
dir=/home/cdaniels/uofc_data/ubs_seq/UBS-seq_basic/data/

cd $dir
rm -f *.gz
ls *.fastq | parallel "seqkit sample -p $prct {} | gzip > test{= s/SRR2353829(.*)\.fastq/\1/; =}_R1.fq.gz"
