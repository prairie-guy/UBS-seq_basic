#!/usr/bin/env python3

# `utils.py`
#
# C. Bryan Daniels, cdaniels@nandor.net
# 03/17/2024

# These functions are generally designed to work with existing `files` and `paths` (in contrast to working
# with sample names generated from`config.yaml` . See `configure.py)
#
# The primary constructors are `mkpath(path)` and 'fname(path,stem,suffix)`
#
# An example of a  general workflow in the pipline might be:
#   - [pwd() -> workspace]
#   - in_path, outpath  = mkpath('path_1'), makepath('path_2'): Create dir withing [workspace]
#   - for sample in samples():
#       ! <bash command> -i {fnames(in_path,sample,'sam')} -o {fname(out_path, sample, 'bam')}
#   altenatively:
#   - for sample in fnames(in_path, 'sam'): # By default only the stem is returned
#       ! <bash command> -i {fname(in_path,sample,'sam')} -o {fname(out_path, sample, 'bam')}
#   or:
#   - for
#   - in_path, outpath  = mkpath('path_2'), makepath('path_3'): Create dir withing [workspace]
#
#   The key difference between samples() and fnames() is that samples() generates its samples
#   from `config.yaml` and fnames() generates its samples from files in an existing path
#
# Reminder: The term `fname` refers to `files` and `path`

import os, sys, re, subprocess, itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from subprocess import getoutput as run
from pathlib import Path
from collections.abc import Iterator
from importlib import reload
from operator import itemgetter as items

def fname(path, stem, suffix, tail=''):
    'Return a `path` and `suffix` complete filename. Optionally, append _`tail` to stem'
    return Path(path)/f"{stem}.{suffix}" if not tail else Path(path)/f"{stem}_{tail}.{suffix}"

def mkpath(path):
    'Return `path` name; creates (not over-writting) a new dir within the pwd. Also prints date/time executed'
    path = Path(path)
    if not os.path.exists(path): os.makedirs(path)
    date = run('date "+%H:%M:%S_%m-%d-%Y"')
    print(f">>> {{{path}}} {date}")
    return path


def fnames(path_or_files, suffix="", key='full'):
    """

    path_or_files(path_or_files, ):: str|Path|list ->
    Returns all files in `path_or_files` with matching `suffix`
    `path_or_files can be either a `path` or a list of `files`
    `stem` is bool. Default: stem=True returns only the stem; stem=False, returns the full filename
    Note: Unlike `fname()`, which creates filenames, `fnames()` returns exisiting files
    """
    if isinstance(path_or_files, (str,Path)):
        if not suffix: suffix = '*'
        files = Path(path_or_files).glob(f'*.{suffix}')
        files = list(map(str,files))

    if isinstance(path_or_files, list):
        if not suffix: suffix = '*'
        files = list(map(str,path_or_files))
        files = list(filter(lambda f: re.search(f'\.{suffix}$',f), files)) # check for extension

    if isinstance(path_or_files, Iterator):
        return fnames(list(path_or_files), suffix, key) # Remember to update signatue

    return list(map(lambda f: fname_index(f)[key], files))


def fnames_string(path_or_files, suffix="", delimiter=' '):
    """
    Returns a string delimited (by delimiter) of files in `path` or in `files` matchiing suffix
    Note: When input are files, no `suffix` is required
    """
    path_or_files =  list(map(lambda f: str(f), path_or_files))
    return  (delimiter.join(path_or_files) if isinstance(path_or_files, list) else
             delimiter.join(list(fnames(path_or_files, suffix))))

def fname_index(path_or_file):
    """
    Returns a dict of a filename split into  keys={'path','keys','name','stem','suffix'}
    """
    path = Path(path_or_file)
    base, _, extension = path.name.partition('.')
    parts = base.split('_')
    return {
        'full': str(path),
        'path': str(path.parent) if path.parent != Path('.') else False,
        'name': path.name,
        'stem': '_'.join(parts),
        'suffix': extension if extension else False,
        'keys': parts if parts else False}

def fname_stem(path_or_file, idx=[0]):
    keys = fname_index(path_or_file)['keys']
    return '_'.join((items(*idx))(keys))


def fnames_match(files, pattern, **wildcards):
    """
    Given a list of `files` and a `pattern` returns a list of list of files matching
    the pattern.

    `pattern` is a string similiar to the target `files` with each wildcard indicated by `{ }`
    `pattern` will generally hava a form: {sample}_{run}_{pair}_{ref}.{ext}, though order
    does not matter. Names do not matter provided they correspond to the names used in
    wildcard named lists.

    `**wilcards` are named lists, where the wildcard name is used within the pattern
    For example: sample = ['t1','c1']

    When files are matched to a pattern, chars between two {wildcard} will be ignored.
    This allows grouping of matched files, which are returned as a list

    files =['t1_1_genes.bam', 't1_r2_genes.bam', 't1_r1_genome.bam', 't1_r2_genome.bam',
            'c1_r1_genes.bam', 'c1_r1_genome.bam']

    match_files(files, "{sample}_{ref}.bam", sample = ['t1','c1'], ref = ['genes','genome']) ->

        [['t1_1_genes.bam', 't1_r2_genes.bam'],
         ['t1_r1_genome.bam', 't1_r2_genome.bam'],
         ['c1_r1_genes.bam'],
         ['c1_r1_genome.bam']]
    """
    def make_regx(pattern):
        "Return a regex pattern, given a pattern with wilds. Chars between {wild}s replaced by .*?"
        pattern = re.sub('}.*?{', r'}.*{', pattern)
        pattern = re.sub('\\.\\w+$', r'\\.\\w+', pattern)
        regex   = re.sub('\\+({\\w+})', r'+\\\\.\\1', pattern)
        return fr"{regex}"

    def expand_wildcards(regx, wildcards):
        "Expands paramaterized `regex` to concrete regular expressions"
        wildcard_names = re.findall(r'{(\w+)}', regx)
        wildcard_values = [wildcards[name] for name in wildcard_names]
        combinations = itertools.product(*wildcard_values)
        expanded_regxes = [regx.format(**dict(zip(wildcard_names, comb))) for comb in combinations]
        #print(expanded_regxes)
        return expanded_regxes

    if isinstance(files,Iterator): files = list(files)
    results = []
    regx = make_regx(pattern)
    #print(regx)
    #print(wildcards)
    for pat in expand_wildcards(regx, wildcards):
        #fn = [m.group(0) for f in files if (m := re.match(pat, f))]
        fn = [f for f in files if re.match(pat, f)]
        if fn: results.append(fn)
    return results

# def make_regx(pattern):
#     pattern = re.sub('}.*?{', r'}.*{', pattern)
#     pattern = re.sub('\.(\w+)$', r'\.\1', pattern)
#     refex   = re.sub('\+(\{ext\})', r'+\\.\1', pattern)
#     return fr"{regex}"

# def make_regx(pattern):
#     pattern = re.sub('}.*?{', r'}.*{', pattern)
#     pattern = re.sub('\\.\\w+$', r'\\.\\w+', pattern)
#     pattern = re.sub('({\\w+})/({\\w+})', r'\\1/\\2', pattern)
#     regex   = re.sub('\\+({\\w+})', r'+\\\\.\\1', pattern)
#     return fr"{regex}"

#
# Utilities to be used with Jupyter
#

def nlines(file):
    'Returns fast linecout (fast)'
    result = subprocess.run(['wc', '-l', file], stdout=subprocess.PIPE)
    n = int(result.stdout.split()[0])
    return n

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

def make_histogram(ds, ds_name, table_label=None, y_label="Frequency", density=True):
    'Makes a histogram for a dataset and its name'
    fig, ax = plt.subplots()
    ax.hist(ds, density=density)
    ax.set_ylabel(y_label)
    ax.set_xlabel(ds_name)
    if table_label is not None:
        ax.set_title(table_label)
    plt.show()
