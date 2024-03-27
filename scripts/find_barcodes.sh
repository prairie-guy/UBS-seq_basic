#!/usr/bin/env sh

# Finds 10 most common 6-mer barcodes immediately preceding an adapter

adapter=$1
file=$2

zcat -c "$file" | grep -o "[ATCG]\{6\}$adapter" | grep -o "^[ATCG]\{6\}" | sort | uniq -c | sort -rn 2> /dev/null| cat 2> /dev/null | head

echo `zcat -c "$file" | grep -o "$adapter" | wc -l` Total Adapters

echo `zcat -c "$file" | wc -l` Total Reads
