#!/usr/bin/env sh

if [ -z "$1" ]; then
  prct=0.01
else
  prct="$1"
fi

cd $(dirname "$(readlink -f "$0")")
cd ../data
rm *.gz 2>/dev/null
ls *.fastq | sort | parallel "seqkit sample -p $prct {} | gzip > test{#}_R1.fq.gz"
