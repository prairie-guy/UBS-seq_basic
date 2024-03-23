#!/usr/bin/env sh

# Changing Random Seed from default of 11;
# Note: As 'random' is deterministic, PE reads will have corresponding pairs selected
seed="131"

if [ -z "$1" ]; then
  prct=0.01
else
  prct="$1"
fi

cd $(dirname "$(readlink -f "$0")")
cd ../data
rm *.gz 2>/dev/null
ls *.fastq | parallel "seqkit sample -p $prct {} | gzip > {.}.fq.gz"
