#!/usr/bin/env python3

# Key Utilities Required By Pipeline
#
#

import os, sys, subprocess, re,itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from importlib import reload

def fname(path, base, suffix, tail=''):
    'Return a `path` and `suffix` complete filename. Optionally, append _`tail` to base'
    return Path(path)/f"{base}.{suffix}" if not tail else Path(path)/f"{base}_{tail}.{suffix}"


# def fname_split(filepath):
#     filepath = Path(filepath)
#     base, _, extension = filepath.name.partition('.')
#     parts = base.split('_')
#     return {
#         'path': filepath.parent if filepath.parent != Path('.') else False,
#         'keys': parts if parts else [],
#         'ext': extension if extension else False}

def fname_split(filepath):
    filepath = Path(filepath)
    base, _, extension = filepath.name.partition('.')
    parts = base.split('_')
    return {
        'path': filepath.parent if filepath.parent != Path('.') else False,
        'keys': parts if parts else [],
        'ext': extension if extension else False,
        'stem': '_'.join(parts)}


def fnames_index(path, suffix):
    """
    Takes a `path` and `suffix` and returns a list of keys delmited by '_' by each fname
    `tail=False` removes the tail key

    fnames_index(in_path,'bam') ->
        [['t1', 'r2', 'genome'],
         ['t1', 'r2', 'genes'],
         ['t1', 'r1', 'genome'],
         ['c1', 'r1', 'genes'],
         ['c1', 'r1', 'genome'],
         ['t1', 'r1', 'genes']]
    """
    pattern = f'*.{suffix}'
    return(p.stem.split('_') for p in Path.glob(path,pattern))


def make_pattern(pattern):
    pattern = re.sub('}.*?{', r'}.*{', pattern)
    pattern = re.sub('\.(\w+)$', r'\.\1', pattern)
    pattern = re.sub('\+(\{ext\})', r'+\\.\1', pattern)
    return fr"{pattern}"


def expand_wildcards(pattern, wildcards):
    wildcard_names = re.findall(r'{(\w+)}', pattern)
    wildcard_values = [wildcards[name] for name in wildcard_names]
    combinations = itertools.product(*wildcard_values)
    expanded_strings = [pattern.format(**dict(zip(wildcard_names, comb))) for comb in combinations]
    return expanded_strings

def match_files(files, pattern, **wildcards):
    results = []
    regx = make_pattern(pattern)
    #print(regx)
    #print(wildcards)
    for pat in expand_wildcards(regx, wildcards):
        fn = [m.group(0) for f in files if (m := re.match(pat, f))]
        if fn: results.append(fn)
    return results

def mkpath(path):
    'Return `path` name; creates (not over-writting) a new dir within the pwd. Also prints date/time executed'
    path = Path(path)
    if not os.path.exists(path): os.makedirs(path)
    date = subprocess.getoutput('date "+%H:%M:%S_%m-%d-%Y"')
    print(f">>> {{{path}}} {date}")
    return path

def nlines(file):
    'Returns fast linecout (fast)'
    result = subprocess.run(['wc', '-l', file], stdout=subprocess.PIPE)
    n = int(result.stdout.split()[0])
    return n

def samples_string(samples,path,suffix='bam'):
    'Returns a space delimited string of sample files'
    return " ".join([str(fname(path,sample,suffix)) for sample in samples])

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
