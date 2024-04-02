# Basic Jupyter Pipeline Template

This is an example of a Jupyter notebook template meant to be edited to create new bioinformatics pipelines. The idea is to make it easy to interactively develop a pipeline, explore the results from each step and then either making adjustments to the input or else moving to the next step. Moreover, the use of graphs and tables at each step can be utilized to better understand both the tools and the data before moving to the next step. Once the pipeline is working as expected, it can be converted to a python script and edited to remove non-essential elements.

Basic assumptions are made about the form of the file structure so that only notebooks, scripts and other elements of the pipeline are saved to github. The large genomic reference files and intermediate workspace directories are ignored by github. If necessary, these files can be saved locally. 

### UBSseq Pipeline Jupyter Notebook
[workspace/ubs_basic.ipynb](UBSseq Pipeline)

### Github File Structure

```
.
├── config.yaml
├── data
│   ├── 41587_2023_2034_MOESM4_ESM.xlsx
│   ├── README_data.md
├── dir_tree
├── environment.yaml
├── img
│   └── adapter.png
├── README.md
├── reference
│   ├── genome
│   ├── index
│   └── README_reference.md
├── samples_info.txt
├── samples.org
├── scripts
│   ├── configure.py
│   ├── downsample.sh
│   ├── find_barcodes.sh
│   ├── fnames.py
│   ├── generate_reference.sh
│   ├── map_se -> ../workspace/map_se
│   ├── nb2py.sh
│   ├── __pycache__
│   ├── README_scripts.md
│   ├── run_ubs_seq.sh
│   └── utils.py
├── Snakefile
├── src
│   └── README_IGV_WEB.md
└── workspace
    ├── calc_methy
    ├── call_converted
    ├── call_filtered_converted
    ├── dedup
    ├── fastqc_post
    ├── fastqc_pre
    ├── filter_calls
    ├── join_pe
    ├── map_se
    ├── merge_se_runs
    ├── README_workspace.md
    ├── select_samples
    ├── trim
    └── ubs_basic.ipynb
```


### .gitignore file

```
reference/
!reference/README_reference.md
data/*
!data/README_data.md
!data/*.py
!data/*.ipynb
workspace/*
!workspace/README_workspace.md
!workspace/*.py
!workspace/*.ipynb
.*~
.snakemake/
```

### The logic for a pipeline is defined through a series of Steps using dirs to save intermediate results. A general workflow might be as follows:

1. Fork this repository from Github and edit as appropriate. (This is a sample for a minimally viable pipeline for UBS-seq.) 
2. For each **Step** in the pipeline, a new dir will be created and labeled **Step** and will contain all files created by that **Step**
3. Within a **Step**, **in_path** and **out_path** will generically refer to the prior and current **Step**
4. Within each **Step**, the required computation occurs. Generally this involves processing files withing **in_path** and saving the resulting files to **out_path**
5. Original samples should be saved in `data/` 
6. **Short, generic filenames** should be symbolically linked to `workspace`/**in_path** from `data/`. The base **filename** will generally be preserved through each **Step** of the pipeline. The dir names and file suffixes will change through each **Step**
7. The function **mkpath(step)** returns a path for a dir **Step**. It will create a dir if need be, but not overwrite an existing dir
8. The function **fname(path,sample,suffix)** returns a file name without actually creating the file


# UBS-seq_basic

The purpose of this project is run a minimally viable UBS-seq pipline. For simplicity, it will run several single-end samples, mapping only to the genome. The core steps of the pipeline are:

- select_samples
- fastqc_pre
- join_pe
- trim
- fastqc_post
- map_se
- merge_se_runs
- dedup
- call_converted
- call_filtered_converted
- calc_methy
- analysis_and_annotation

This pipeline is based upon the paper by [Qing Dai, etal](https://doi.org/10.1038/s41587-023-02034-w) and the UBS-seq pipeline developed by [Chang Ye](https://github.com/y9c/m5C-UBSseq)


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
- cd data
- ls *.fastq | parallel 'seqkit sample -p 0.01 {} | gzip > test{= s/SRR2353829(.*)\.fastq/\1/; =}_R1.fq.gz'


### workflow
