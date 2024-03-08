#!/usr/bin/env bash

cd $(dirname "$(readlink -f "$0")")
cd ../workspace/
rm -fr call_converted  call_filtered_converted  dedup  fastqc_post  fastqc_pre  filter_calls  map_genes  map_genome samples trim
