#!/usr/bin/env bash

cd $(dirname "$(readlink -f "$0")")
cd ../workspace/
rm -fr select_samples fastqc_pre join_pe trim call_converted call_filtered_converted combined dedup fastqc_post map_pe map_se merge_pe_runs merge_se_runs
