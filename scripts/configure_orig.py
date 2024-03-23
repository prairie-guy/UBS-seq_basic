#!/usr/bin/env python3

# `configure.py`
#
# C. Bryan Daniels, cdaniels@nandor.net
# 03/17/2024
# sample2data and group2sample data structures developed by Chang Ye


"""
`configure.py`

Configure sequencing samples from config.yaml

The primary data structures created from config.yaml are:
    - sample2data

    - group2sample

    - group2run

Major functions and operators operating on these data structures are:
    - samples(n=2, se=True, keys=False, extra=[])

    - data(se=True, keys=False)

    - groups()

These operators are NOT intended to be used for actual files. (Use fnames() in fnames.py instead)

    - Create data structures and functions to operate on samples defined in confi.yaml

    - This module can consume and generate data structures as defined in config.yaml

    - Generally, the return type is iterator

(See `fnames.py` for operators to be used for use with paths and files)

ToDo: Add additional functionality for other elements of `config.yaml`, including `references`
"""


from pathlib import Path
from snakemake import load_configfile
from collections import defaultdict
from itertools import compress

home_path      = Path.cwd()/'..'
config = load_configfile(home_path/"config.yaml")

def parse_samples(config):
    #read_ids = ["R1", "R2"]
    pairend_run_ids = []
    sample2data, group2sample, group2run  = defaultdict(dict), defaultdict(list), defaultdict(list)
    for s, v in config["samples"].items():
        group2sample[v["group"]].append(s)
        for i, v2 in enumerate(v["data"], 1):
            r = f"r{i}"
            group2run[v["group"]].append(f"{s}_{r}")
            sample2data[s][r] = {k3: home_path/v3 for k3, v3 in v2.items()}  # Using string representation for simplicity
            if len(v2) == 2:
                pairend_run_ids.append(s + "_" + r)
    return sample2data, group2sample, group2run, pairend_run_ids


# Defines Sample Data Structures for Export
sample2data, group2sample, group2run, pairend_run_ids = parse_samples(config)

def samples(n=2, se=True, keys=False, extra=[]):
    '''
    samples :: Int -> Bool -> Bool -> [String] ->
    Generates iterators of sample names as defined in 'sample2data'
    - n=1 returns an iterator of biological replicates, sample(2) -> ['t1', 't2', 'c1']
    - n=2 (default) returns an iterator of biologial replicates and runs -> ['t1_r1', 't2_r1', 'c1_r1']
    - se=True (default) filters for SE reads, se=False filters for PE reads
    - keys=True returns an iterator of tuples of length 'n'. When iterated will require 'n' keys -> [('t1', 'r1'), ('t2', 'r1'), ('c1', 'r1')]
    - extra=[alist] returns an iterator which is the Cartesian product of the elements of alist and the original iterator

    Example:
    samples(se=True, keys=True, extra=['gene','genome']) ->
     [('t1', 'r1', 'gene'), ('t1', 'r1', 'genome'), ('t2', 'r1', 'gene'), ('t2', 'r1', 'genome'), ('c1', 'r1', 'gene'), ('c1', 'r1', 'genome')]
    '''
    if n not in [1, 2]: return("usage: samples(n=2, keys=False, extra=False)")
    if not isinstance(se,bool) or not isinstance(keys,bool): return("usage: samples(n=2, keys=False, extra=False)")
    valid_extra = isinstance(extra, list) and len(extra) > 0
    # Generators
    # [n,keys,extra]
    results = {
        1: {(False, False): (s for s in sample2data.keys()),
            (True, False):  (s for s in sample2data.keys()),
            (False, True):  (f"{s}_{e}" for s in sample2data.keys() for e in extra),
            (True, True):   ((s, e) for s in sample2data.keys() for e in extra)},

        2: {(False, False): (f"{s}_{r}" for s, v in sample2data.items() for r in v.keys()),
            (True, False):  ((s, r) for s, v in sample2data.items() for r in v.keys()),
            (False, True):  (f"{s}_{r}_{e}" for s, v in sample2data.items() for r in v.keys() for e in extra),
            (True, True):   ((s, r, e) for s, v in sample2data.items() for r in v.keys() for e in extra)}}
    gen = results[n].get((keys, valid_extra), "usage: samples(n=2, keys=False, extra=False)")
    gen = list(gen);pprint(gen);gen=(gen)
    if n ==1:
        return gen
    # Masks to filter SE and PE Reads
    if n ==2:
        # [se][extra]
        results = {
            True:  {True:   [len(sample2data[s][r])==1 for s in sample2data for r in sample2data[s]],
                    False:  [len(sample2data[s][r])==1 for s in sample2data for r in sample2data[s] for e in extra]},
            False: {True:   [len(sample2data[s][r])==2 for s in sample2data for r in sample2data[s]],
                    False:  [len(sample2data[s][r])==2 for s in sample2data for r in sample2data[s] for e in extra]}}

        mask = results[se][not extra]

        if se: return compress(gen,mask)
        if keys: return ((s,r,'R1','R2') for s,r in compress(gen,mask))
        return ((f'{k}, ''R1','R2') for k in compress(gen,mask))

def data(se=True, keys=False):
    """
    data :: Bool -> Bool -> [(String, Path)] | [(String, Path, Path)]
    Generates iterator of data paths as defined in `sample2data`
    - type='SE' will flatten the list and return a key[s] with a single path to R1 sequence reads
    - type='PE' will flatten the list and return a key[s] with a two paths to R1,R2  sequence reads

    Example:
    data() -> [('t1_r1', PosixPath('../data/test1_R1.fq.gz')),
               ('t2_r1', PosixPath('../data/test2_R1.fq.gz')),
               ('c1_r1', PosixPath('../data/test3_R1.fq.gz'))]

    data(se=False) ->
              [('t1_r1', PosixPath('../data/test1.1_R1.fq.gz'), PosixPath('../data/test1.1_R2.fq.gz')),
               ('t1_r2', PosixPath('../data/test_1.2_R1.fq.gz'), PosixPath('../data/test1.2_R2.fq.gz')),
               ('t1_r3', PosixPath('../data/test1.3_R1.fq.gz'), PosixPath('../data/test1.3_R2.fq.gz'))]
    """
    # if type == 'SE': return ((f"{s}_{r}",d) for s, v in sample2data.items() for r, ds in v.items() for d in ds.values())
    # if type == 'PE': return ((f"{s}_{r}", ds) for s,v in sample2data.items() for r,ds in v.items())
    # return "usage: data(type), type = {'SE' | 'PE'}"
    if not isinstance(se,bool): return("usage: data(se=True), for PE  data(se=False)")
    # [se][keys]
    results = {
        True:  {False: ((f'{s}_{r}',sample2data[s][r]['R1']) for s,r in samples(se=True,keys=True)),
                True:  ((s,r,sample2data[s][r]['R1']) for s,r in samples(se=True,keys=True))},
        False: {False: ((f'{s}_{r}', sample2data[s][r][r1], sample2data[s][r][r2]) for s,r,r1,r2 in samples(se=False,keys=True)),
                True:  ((s,r, sample2data[s][r][r1], sample2data[s][r][r2]) for s,r,r1,r2 in samples(se=False,keys=True))}}
    return results[se][keys]

def groups(n=2):
    if n == 1: return ((s, g) for g, items in group2sample.items() for s in items)
    if n == 2: return ((s, g) for g, items in group2run.items() for s in items)
