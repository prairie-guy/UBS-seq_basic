#!/usr/bin/env bash

cd $(dirname "$(readlink -f "$0")")
cd ../workspace/
rm -fr fastqc_post/ fastqc_pre/ hisat3n_align/ hisat3n_call/ hisat3n_dedup/ hisat3n_sort/ trim/

stamp=`date "+%d-%m-%Y-%H:%M"`

ipython ubs_basic.py | tee $stamp.log
