#!/usr/bin/env sh

pd=`pwd`

genome=$pd/reference/genome
index=$pd/reference/index/hisat3n

gen_cont=$genome/contamination
gen_snc=$genome
gen_gen=$genome

idx_cont=$index/contamination
idx_snc=$index/Homo_sapiens.GRCh38.sncRNA
idx_gen=$index/Homo_sapiens.GRCh38.genome

build_index=$pd/src/hisat-3n/hisat-3n-build

mkdir -p $gen_cont $idx_cont $idx_snc $idx_gen

# Download fasta files
cd $gen_cont
wget https://github.com/chelab/db/raw/main/reference_sequence/contamination.fa.gz
gunzip contamination.fa.gz
cd $pd

cd $gen_snc
wget https://github.com/chelab/db/raw/main/reference_sequence/Homo_sapiens.GRCh38.sncRNA.fa.gz
gunzip Homo_sapiens.GRCh38.sncRNA.fa.gz
cd $pd

cd $gen_gen
wget http://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz Homo_sapiens.GRCh38.genome.fa.gz
gunzip Homo_sapiens.GRCh38.genome.fa.gz
cd $pd

# Generate indexes
cd $idx_cont
$build_index --base-change C,T $gen_cont/contamination.fa contamination
cd $pd

cd $idx_snc
$build_index --base-change C,T $gen_snc/Homo_sapiens.GRCh38.sncRNA.fa Homo_sapiens.GRCh38.sncRNA
cd $pd


cp $pd/Homo_sapiens.GRCh38.genome.fa $gen_gen/. ### REM to REMOVE
cd $idx_gen
$build_index --base-change C,T $gen_gen/Homo_sapiens.GRCh38.genome.fa Homo_sapiens.GRCh38.genome
cd $pd
