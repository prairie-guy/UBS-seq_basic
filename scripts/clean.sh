#!/usr/bin/env bash

cd $(dirname "$(readlink -f "$0")")
cd ../workspace/
rm -fr select_samples trim fastqc_post  fastqc_pre  map_se  map_se_pe map_combined call_converted  call_filtered_converted  dedup
