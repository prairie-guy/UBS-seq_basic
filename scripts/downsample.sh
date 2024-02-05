#!/usr/bin/env sh

if [ -z "$1" ]; then
  prct=0.01
else
  prct="$1"
fi

cd $(dirname "$(readlink -f "$0")")
cd ../data
rm -f *.gz
ls *.fastq | parallel "seqkit sample -p $prct {} | gzip > test{= s/SRR2353829(.*)\.fastq/\1/; =}_R1.fq.gz"
