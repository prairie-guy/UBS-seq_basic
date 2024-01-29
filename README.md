# UBS-seq_basic

### Directories and Files added to .gitignore
- data/
- reference/
- scr/
- .*~ 

### How to create `reference/`

#### From `config.yaml`

reference:
  contamination:
    fa: ~/reference/genome/contamination/contamination.fa
    hisat3n: ~/reference/index/hisat3n/contamination/contamination
  genes:
    fa: ~/reference/genome/Homo_sapiens.GRCh38.sncRNA.fa
    hisat3n: ~/reference/index/hisat3n/Homo_sapiens.GRCh38.sncRNA/Homo_sapiens.GRCh38.sncRNA
  genome:
    fa: ~/reference/genome/Homo_sapiens.GRCh38.genome.fa
    hisat3n: ~/reference/index/hisat3n/Homo_sapiens.GRCh38.genome/Homo_sapiens.GRCh38.genome

#### Install histat-3n from source 
- mkdir src
- cd src
- git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n
- cd hisat-3n
- git checkout -b hisat-3n origin/hisat-3n
- make
- binaries at: src/histat-3n/

#### Run script `generate_reference.sh`
- Assumes histat-3n-build located at:  src/histat-3n/histat-3n-build
- rm -fr reference (This will not fix a partially installed reference/)
- Run from main repo
- scripts/generate_reference.sh

### Downsize Samples to ~ 10 milllion reads per file
- ls *.fastq | parallel 'seqkit sample -p 0.01 {} | gzip > test{= s/SRR2353829(.*)\.fastq/\1/; =}_R1.fq.gz'


### workflow

![](./docs/flow.svg)
