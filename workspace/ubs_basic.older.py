#!/usr/bin/env python
# coding: utf-8

# # UBS-seq Pipeline
# ## Basic Workflow

# **C. Bryan Daniels**
# 
# **2/1/2024**

# ## Project: UBS-seq Basic Pipeline

# The purpose of this project is run a minimally viable UBS-seq pipline. For simplicity, it will run several single-end samples, mapping only to the genome. The core steps of the pipeline are:
# - cut_apapter
# - quality_control
# - align2ref
# - sort2ref
# - dedupe
# - filter->all_multi_unique
# - call_peaks
# - select_groups
# - analysis_and_annotation
# 
# This pipeline is based upon the paper by [Qing Dai, etal](https://doi.org/10.1038/s41587-023-02034-w) and the UBS-seq pipeline developed by [Chang Ye](https://github.com/y9c/m5C-UBSseq)
# 
# 

# ## Setup

# #### The logic for the Pipeline is defined through a series of Steps using dirs to save intermediate results
# 1. For each **Step** in the pipeline a dir will be created and labeled **Step** and will contain all files created by that **Step**
# 2. Within a **Step**, **in_path** and **out_path** will generically refer to the prior and current **Step**
# 3. Within each **Step**, the appropriate processes will occur. Generally this involves processing files from **in_path** and saving to **out_path**
# 4. **Abbreviated filenames** should not change through the pipeline (suffixes will reflect current file formats). The dir name should reflect the **Step**, not the filename.
# 6. The function **mkpath(step)** returns a path for a dir **Step**. It will create a dir if need be, but not overwrite an existing dir
# 8. The function **fname(path,sample,suffix)** returns a file name without actually creating the file

# #### Execution from Command Line

# - cd workplace/
# - juptyer nbconvert ubs_basic.ipynb --to script
# - ipython ubs_basic.py

# #### Environment

# In[ ]:


import os, sys, subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from pathlib import Path
from IPython.display import display, HTML
from snakemake import load_configfile


# In[ ]:


def fname(path, base, sufix):
    'Return a path and suffix complete filename'
    return Path(path)/f"{base}.{sufix}"

def mkpath(path):
    'Return dir path name; creates (not over-writting) a new dir within the pwd. Also prints date/time executed'
    path = Path(path)
    if not os.path.exists(path): os.makedirs(path)
    date = get_ipython().getoutput(' date "+%H:%M:%S_%m-%d-%Y"')
    print(f">>> {{{path}}} {date[0]}")
    return path


# #### Path References

# `Waring:` This notebook should be located and executed from within the directory `home_path/workspace/`
# 
# `Note:` Github ignores the data/, workspace/ and reference directories, except for `.py` and `.ipynb` files

# In[ ]:


home_path      = Path.cwd()/'..'
data_path      = home_path/'data'
workspace_path = home_path/'workspace'


# Use `config.yaml` to configure `references/`, but not samples in `data/`

# In[ ]:


config = load_configfile("../config.yaml")


# In[ ]:


genome_fa  = home_path/config['reference']['genome']['fa'].removeprefix('~/')
genome_idx = home_path/config['reference']['genome']['hisat3n'].removeprefix('~/')


# In[ ]:


#genome_fa = ref_path/'genome/Homo_sapiens.GRCh38.genome.fa'
#genome_idx = ref_path/'index/hist3n/Homo_sapiens.GRCh38.genome'


# In[ ]:


# Add to shell PATH
os.environ['PATH'] = f"{str(home_path)}:" + os.environ['PATH'] # home_path
os.environ['PATH'] = '/home/cdaniels/bin/homer:' + os.environ['PATH'] # homer
os.environ['PATH'] = '/home/cdaniels/bin/hisat-3n:' + os.environ['PATH'] # hisat-3n


# In[ ]:


# Set Java Flag
os.environ['_JAVA_OPTIONS'] = '-Xmx8g'


# In[ ]:


# Number of cores                                                                                                                                                                                               
nc = get_ipython().getoutput('nproc')                                                                                                                                                                           
nc = int(nc[0])                                                                                                                                                                                                 
nc  


# #### Functions

# In[ ]:


def nlines(file):
    'Returns fast linecout (fast)'
    result = subprocess.run(['wc', '-l', file], stdout=subprocess.PIPE)
    n = int(result.stdout.split()[0])
    return n


# In[ ]:


def nseqs(bam_fastq):
    'Returns number of sequences in bam, sam, fasta or fastq file'
    n = get_ipython().getoutput('samtools view -c {bam_fastq}')
    return int(n[0])


# In[ ]:


def samples_string(samples,path,suffix='bam'):
    'Returns a space delimited string of sample files'
    return " ".join([str(fname(path,sample,suffix)) for sample in samples])    


# In[ ]:


def make_table(ds1, ds2, ds1_name, ds2_name, y_label=None, xs_labels=None, table_label=None):
    'Makes a bar graph comparing two datasets and their corresponding names'
    assert(len(ds1)==len(ds2))
    n = np.arange(len(ds1))
    width = 0.35
    fig, ax = plt.subplots()
    rects1 = ax.bar(n - width/2, ds1, width, label=ds1_name)
    rects2 = ax.bar(n + width/2, ds2, width,label=ds2_name)
    ax.set_ylabel(y_label)
    ax.set_title(table_label)
    ax.set_xticks(n)
    ax.set_xticklabels(range(1,len(ds1)+1)) if xs_labels == None else ax.set_xticklabels(xs_labels)
    ax.legend()
    plt.show()


# In[ ]:


def make_histogram(ds, ds_name, table_label=None, y_label="Frequency", density=True):
    'Makes a histogram for a dataset and its name'
    fig, ax = plt.subplots()
    ax.hist(ds, density=density)
    ax.set_ylabel(y_label)
    ax.set_xlabel(ds_name)
    if table_label is not None:
        ax.set_title(table_label)
    plt.show()


# ## Step: select_samples

# #### Select Samples and Filenames

# In[ ]:


out_path = mkpath('samples')


# #### Project definitions for treated samples, control samples

# In[ ]:


# short sample names
treated = ["t1","t2","t3"] # Edit
control = [] # Edit
samples = treated + control


# In[ ]:


print(samples)


# In[ ]:


# sample filenames
samples_fn = [data_path/'test0_R1.fq.gz', data_path/'test1_R1.fq.gz', data_path/'test2_R1.fq.gz']


# In[ ]:


# dict
s2fn = {name: fname for (name,fname) in zip(samples,samples_fn)};
s2fn


# In[ ]:


for (s,fn) in s2fn.items():
    get_ipython().system(' ln -s {data_path/fn} {fname(out_path,s,"fq.gz")}')


# Check that the files look correct

# In[ ]:


get_ipython().system('ls -lLh {out_path}')


# How many sequence reads do we have per Sample

# In[ ]:


# Total Reads per Samples

#files = [fname(out_path,sample, "fq.gz") for sample in samples]
#for f in files:
#    ns = nseqs(f)
#    print(f"{f}: {ns:,}")


# ## Step: fastqc_pre

# #### Pre Trimming Quality Control

# In[ ]:


in_path = mkpath("samples")
out_path = mkpath("fastqc_pre")


# #### fastqc

# In[ ]:


#! fastqc --help


# ##### paramaters:
# - -o output dir

# In[ ]:


for sample in samples:
    get_ipython().system(' fastqc {fname(in_path,sample,"fq.gz")} -o {out_path} 2> /dev/null')


# In[ ]:


get_ipython().system(' ls {out_path}')


# ## Step: trim

# #### Trim the adapter and downstream sequence as well as trimmng lower quality downstream sequence

# In[ ]:


in_path = mkpath("samples")
out_path = mkpath("trim")


# In[ ]:


adapter =  "AGATCGGAAGAGCACACGTCT"
barcode3 = "ATCACG"


# In[ ]:


get_ipython().system('ls {in_path}')


# #### cutadapt

# In[ ]:


#! cutadapt --help


# ```
# cutadapt -j {threads} 
#             -n 2 
#             -a "{params.barcode3}{params.adapter3};e=0.15;o=6;anywhere;" 
#             --untrimmed-output={output.fastq_untrimmed} 
#             -o - {input} 2>{output.report1} | 
# cutadapt -j {threads} 
#             -u 5 -u -5 
#             --rename='{{id}}_{{cut_prefix}}{{cut_suffix}} {{comment}}' 
#             --max-n=0 
#             -q 15 
#             --nextseq-trim=15 
#             -m 20 
#             --too-short-output={output.fastq_tooshort} 
#             -o {output.fastq_cut} - >{output.report2}
# ```

# https://cutadapt.readthedocs.io/en/stable/guide.html

# ##### paramaters:
# - -j 0, Number of CPU. Use 0 to autodetect
# - --nextseq-trim=15, is used to trim these low-quality tails of 'G's by introducing a variable quality threshold
# - --action=trim, trim adapter and downstream sequence
# - -a '{barcode3}{adapter};e=0.15;o=6;anywhere;', 
# - -n 2, Remove up to COUNT adapters from each read. Default: 1
# - -u 5 -u -5, Remove LEN bases from each read. If positve, from beginning. If negative, from end
# - --max-n=0, Discard reads with more than COUNT 'N' bases (here any)
# - -q 15, Trim low-quality bases from 5' end of eachread before adapter removal.
# - -m 20, Discard reads shorter than LEN
# - --rename='{{id}}_{{cut_prefix}}{{cut_suffix}} {{comment}}', rename comment of fastq record
# - --too-short-output={fname(out_path,sample,"fastq_tooshort")},
# - -o {fname(out_path,sample,"fq.gz")}, output file
# {fname(in_path,sample,"fq.gz")}, input file

# **NOTE:** Added `--length 30` to cutoff reads where C begin to become more common. Not in original code

# In[ ]:


for sample in samples:
    get_ipython().system('cutadapt -j 0 --nextseq-trim=15 --action=trim -a \'{barcode3}{adapter};e=0.15;o=6;anywhere;\'             -n 2 -u 5 -u -5 --max-n=0 -q 15 -m 20 -l 80             --length 30              --rename=\'{{id}}_{{cut_prefix}}{{cut_suffix}} {{comment}}\'             --too-short-output={fname(out_path,sample,"fastq_tooshort")}              -o {fname(out_path,sample,"fq.gz")}               {fname(in_path,sample,"fq.gz")} > {fname(out_path,sample,"log")}')


# #### Analysis

# Take a look at a FASTQ file

# In[ ]:


in_fn  = fname(in_path,treated[0],'fq.gz')
out_fn = fname(out_path,treated[0],'fq.gz')
out_fn


# In[ ]:


get_ipython().system(' zcat {in_fn}|head -16')


# In[ ]:


def show_adapter(reads):
    for read in reads.split('\n'):
        read = read.replace(adapter, f'<span style="color: blue;">{barcode3}{adapter}</span>')
        display(HTML(read))


# Look for adapter in untrimmed reads

# In[ ]:


reads = get_ipython().getoutput(" zcat {in_fn} | head -36  | seqtk seq -A |grep -v '>'")
reads = ('\n').join(reads)
show_adapter(reads)


# Verify that adapters and all downstream elements of reads have been trimmed

# In[ ]:


reads = get_ipython().getoutput("zcat {out_fn}| head -36  | seqtk seq -A |grep -v '>'")
reads = ('\n').join(reads)
show_adapter(reads)


# In[ ]:


get_ipython().system('ls {out_path}')


# What did `--rename='{{id}}_{{cut_prefix}}{{cut_suffix}} {{comment}}'` do?
# 
# It appears to has included in the fastq header, an id comprised of the first and last 5 nt cut by `-u 5` and `-u -5`

# In[ ]:


get_ipython().system(' zcat {in_fn}|head -2')


# In[ ]:


get_ipython().system('zcat {out_fn}| head -2')


# Too slow for large files. Commented out.

# In[ ]:


#files = [fname(in_path,sample, "fq.gz") for sample in samples]
#res = []
#for f in files:
#    n = !seqtk seq -A {f}|grep -v '>'|wc -c
#    res.append(int(n[0]))
#ins = res
#ins


# In[ ]:


#files = [fname(out_path,sample, "fq.gz") for sample in samples]
#res = []
#for f in files:
#    n = !seqtk seq -A {f}|grep -v '>'|wc -c
#    res.append(int(n[0]))
#outs = res
#outs


# In[ ]:


#make_table(ins, outs, "Origs", "Trimmed", "Sum of Total Read Lengths", samples, "Sum of Total Reads Per Sample")


# In[ ]:


#reads = !xargs zcat {out_path}/*.gz | seqtk seq -A  |grep -v ">" 
#ds = [len(read) for read in reads]
#make_histogram(ds, "Read Length", "Trimmed Reads")


# ## Step: fastqc_post

# #### Post Trimming Quality Control

# In[ ]:


in_path = mkpath("trim")
out_path = mkpath("fastqc_post")


# In[ ]:


for sample in samples:
    get_ipython().system(' fastqc {fname(in_path,sample,"fq.gz")} -o {out_path} 2> /dev/null')


# #### Consolidate fastqc Reports

# In[ ]:


get_ipython().system('multiqc -f -fp -m fastqc -n multiqc -o {out_path} {out_path}')


# In[ ]:


get_ipython().system(' ls {out_path}')


# ## Step: Hisat3n_align

# #### Align Samples to Genome with Hisat-3n

# In[ ]:


in_path = mkpath("trim")
out_path = mkpath("hisat3n_align")


# #### hisat-3n

# http://daehwankimlab.github.io/hisat2/hisat-3n/

# In[ ]:


#! hisat-3n --help


# ```
# hisat3n --index {params.index}
#         -p {threads}
#         --summary-file {output.summary}
#         --new-summary
#         -q
#         -U {input}
#         --directional-mapping
#         --base-change C,T
#         --pen-noncansplice 20
#         --mp 4,1
#         --un {output.fq}
#         -S {output.sam}
# ```

# ##### paramaters:
# - --index {params.index}, Index filename prefix (minus trailing .X.ht2) 
# - -p {nc}, number of threads
# - --summary-file, print alignment summary to this file.
# - --new-summary, print alignment summary in a new style, which is more machine-friendly.
# - -q, query input files are FASTQ .fq/.fastq (default)
# - -U {input}, 
# - --directional-mapping, make directional mapping, please use this option only if your reads are prepared with a strand specific library (off)
# - --base-change C,T,  the converted nucleotide and converted to nucleotide (C,T)
# - --pen-noncansplice 20, penalty for a non-canonical splice site (12)
# - --mp 4,1,  max and min penalties for mismatch; lower qual = lower penalty <6,2>
# - --un {output.fq}, write unpaired reads that didn't align to <path> 
# - -S {output.sam}, File for SAM output (default: stdout)

# In[ ]:


for sample in samples:
    get_ipython().system('hisat-3n --index {genome_idx}         -p {nc}         --summary-file {fname(out_path,sample,"summary")}         --new-summary         -q         -U {fname(in_path,sample,\'fq.gz\')}         --directional-mapping         --base-change C,T         --pen-noncansplice 20         --mp 4,1         --un {fname(out_path,sample,\'unmapped.fq\')}         -S {fname(out_path,sample,\'sam\')}')


# #### Analysis

# In[ ]:


get_ipython().system('ls -lh {out_path}')


# In[ ]:


get_ipython().system(' cat {out_path}/t1.summary')


# In[ ]:


get_ipython().system(' head -8 {out_path}/t1.unmapped.fq')


# In[ ]:


get_ipython().system(" grep -v '@' {out_path}/t1.sam |head -1")


# ## Step: Hisat3n_sort

# #### Sort and Index Hisat3n Sam Files

# In[ ]:


in_path = mkpath("hisat3n_align")
out_path = mkpath("hisat3n_sort")


# #### samtools

# http://www.htslib.org/doc/samtools.html

# ```
# samtools view
#     -@ {threads}
#     -F4 -b {input} |
# samtools sort
#     -@ {threads}
#     --write-index
#     -m 4G
#     -O BAM
#     -o {output} -
# ```

# ##### paramaters:
# ```
# samtools view, SAM<->BAM<->CRAM conversion
#     -@ {nc}, number of threads 
#     -F4, have none of the FLAGs present (-F 4 filters out unmapped reads) 
#     -b, output a bam file
#     {input}  |
# samtools sort, sort alignment file
#     -@ {nc}, number of threads
#     --write-index, index the output files
#     -O BAM, output file format
#     -o {output} -
# ```

# In[ ]:


for sample in samples:
    get_ipython().system(" samtools view -@ {nc} -F4 -b {fname(in_path, sample,'sam')}  |        samtools sort -@ {nc} --write-index -O BAM -o {fname(out_path,sample,'bam')} -")


# #### Analysis

# In[ ]:


get_ipython().system(' ls -lh {out_path}')


# ## Step: Hisat3n_dedup

# #### Remove Dulpicate Reads

# In[ ]:


in_path = mkpath("hisat3n_sort")
out_path = mkpath("hisat3n_dedup")


# #### umicollapse

# https://github.com/Daniel-Liu-c0deb0t/UMICollapse

# ##### paramaters:
# ```
# umicollapse bam, use a bam file
#     --two-pass, use a separate two-pass algorithm for SAM/BAM deduplication.
#     -i {input.bam}, indexed input bam file
#     -o {output.bam}, output bam file
#       > {output.log}
# ```

# In[ ]:


for sample in samples:
    get_ipython().system("umicollapse bam       --two-pass      -i {fname(in_path,sample,'bam')}      -o {fname(out_path,sample,'bam')}     >  {fname(out_path,sample,'log')}")


# #### Analysis

# In[ ]:


get_ipython().system(' ls -lh {out_path}')


# In[ ]:


get_ipython().system(' cat {out_path}/t2.log')


# ## Step: Hisat3n_call

# #### Call Converted bases

# In[ ]:


in_path = mkpath("hisat3n_dedup")
out_path = mkpath("hisat3n_call")


# #### hisat-3n-table

# http://daehwankimlab.github.io/hisat2/hisat-3n/

# ##### paramaters:
# ```
# samtools view -e "rlen<100000" -h {input} |
# hisat3ntable
#     -p {threads}
#     -m --alignments -
#     --ref {params.fa}
#     --output-name /dev/stdout
#     --base-change C,T                      |
#     bgzip -@ {threads} -c > {output}
# ```

# In[ ]:


for sample in samples:
    get_ipython().system('samtools view -e "rlen<100000" -h {fname(in_path,sample,\'bam\')} |     hisat-3n-table         -p {nc}         -m --alignments -         --ref {genome_fa}         --output-name /dev/stdout         --base-change C,T                          |     bgzip          -@ {nc}          -c > {fname(out_path,sample,\'tsv.gz\')}')


# #### Analysis

# In[ ]:


get_ipython().system(' ls -lh {out_path}')


# There are 7 columns in the 3N-conversion-table:
# 
# 1. `ref:` the chromosome name.
# 2. `pos:` 1-based position in ref.
# 3. `strand:` ‘+’ for forward strand. ‘-‘ for reverse strand.
# 4. `convertedBaseQualities:` the qualities for converted base in read-level measurement. Length of this string is equal to the number of converted Base in read-level measurement.
# 5. `convertedBaseCount:` number of distinct read positions where converted base in read-level measurements were found. this number should equal to the length of convertedBaseQualities.
# 6. `unconvertedBaseQualities:` the qualities for unconverted base in read-level measurement. Length of this string is equal to the number of unconverted Base in read-level measurement.
# 7. `unconvertedBaseCount:` number of distinct read positions where unconverted base in read-level measurements were found. this number should equal to the length of unconvertedBaseQualities.

# In[ ]:


get_ipython().system(' zcat {out_path}/t1.tsv.gz |head -20')


# In[ ]:


for sample in samples:
    df = pd.read_csv(fname(out_path, sample,'tsv.gz'), sep='\t', compression='gzip', low_memory=False, dtype={4: int, 6: int })
    df_f = df[df['unconvertedBaseCount'] > df['convertedBaseCount']]
    df_f.to_csv(fname(out_path,sample,'called.csv'),index=False)  # Set index=False i

