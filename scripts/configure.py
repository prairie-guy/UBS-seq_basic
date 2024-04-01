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

import re
from pathlib import Path
from snakemake import load_configfile
from snakemake.io import expand
from collections import defaultdict
from cytoolz import concat, unique

home_path      = Path.cwd()/'..'
config = load_configfile(home_path/"config.yaml")

def parse_samples(config):
    pairend_run_ids = []
    sample2data, group2sample, group2run, group2data  = defaultdict(dict), defaultdict(list), defaultdict(list),defaultdict(list)
    for s, v in config["samples"].items():
        group2sample[v["group"]].append(s)
        for i, v2 in enumerate(v["data"], 1):
            r = f"r{i}"
            group2run[v["group"]].append(f"{s}_{r}")
            sample2data[s][r] = {k3: home_path/v3 for k3, v3 in v2.items()}
            k3 = list(v2.keys())
            if len(k3) ==1:
                group2data[v["group"]].append(f"{s}_{r}")
            if len(k3) == 2:
                group2data[v["group"]].append(f"{s}_{r}_{k3[0]}")
                group2data[v["group"]].append(f"{s}_{r}_{k3[1]}")
                pairend_run_ids.append(s + "_" + r)

    sample2list = [(sk,rk,list(rv.values())) for sk,sv in sample2data.items() for rk,rv in sv.items()]
    return sample2list, sample2data, group2sample, group2run, group2data, pairend_run_ids


# Defines Sample Data Structures for Export
sample2list, sample2data, group2sample, group2run, group2data, pairend_run_ids = parse_samples(config)


def samples(n=2, end='runs', extra=[]):
    """
    samples :: Int -> Bool -> [String] ->
    Generates iterators of sample names as defined in 'sample2data'
    - n=1 returns biological replicates, sample(2) -> ['t1', 't2', 'c1']
    - n=2 (default) biologial replicates with run -> 't1_r1', 't1_r2', 't2_r1', 't2_r2', 'c1_r1']
    - end='se' filters for SE reads -> ['t1_r1', 't1_r2', 'c1_r1']
    - end='pe' filters for pairs of PE reads -> [('t2_r1_R1', 't2_r1_R2'), ('t2_r2_R1', 't2_r2_R2')]
    - end='all' returns all (flat) -> ['t1_r1', 't1_r2', 'c1_r1', 't2_r1_R1', 't2_r1_R2', 't2_r2_R1', 't2_r2_R2']
    - end='runs' biologial replicates by run -> (i.e., runs) -> ['t1_r1', 't1_r2', 't2_r1', 't2_r2', 'c1_r1']
    - end='pe_runs' pe biological replicates by run -> ['t2_r1', 't2_r2']
    - extra=[alist] returns an iterator which is the Cartesian product of the elements of alist and the original iterator

    Example:
    samples(end='se',extra=['gene','genome']) ->
     [('t1', 'r1', 'gene'), ('t1', 'r1', 'genome'), ('t2', 'r1', 'gene'), ('t2', 'r1', 'genome'), ('c1', 'r1', 'gene'), ('c1', 'r1', 'genome')]
    """
    if n not in [1, 2]: return("usage: samples(n=2, end='se', extra=[])")
    if end not in ['se','pe','all','runs','pe_runs'] : return("usage: samples(n=2, end='se', extra=[])")
    if not isinstance(extra, (list,tuple)): return("usage: samples(n=2, end='se', extra=[])")

    # Same for end = 'se'|'pe'
    if n == 1:
        result = list(sample2data.keys())
        return result
        return result if not extra else [f'{s}_{e}' for s in result for e in extra]

    if end == 'se':
        result = [f'{s}_{r}' for s,r,d in sample2list if len(d)==1]
        return result if not extra else [f'{s}_{e}' for s in result for e in extra]

    if end == 'pe':
        result = [(f'{s}_{r}_R1',f'{s}_{r}_R2') for s,r,d in sample2list if len(d)==2]
        return result if not extra else [(f'{s1}_{e}',f'{s2}_{e}') for s1,s2 in result for e in extra]

    if end == 'all':
        return samples(n=n,end='se', extra=extra) + list(concat(samples(n=n,end='pe', extra=extra)))

    if end == 'runs':
        result = [f'{s}_{r}' for s,r,d in sample2list]
        return result if not extra else [f'{s}_{e}' for s in result for e in extra]

    if end == 'pe_runs':
        result = list(concat(samples(end='pe',extra=extra)))
        result =list(map(lambda string: re.sub(r'_(?:R1|R2)(?=_|$)', '', string),
                        result))
        return list(unique(result))

    return("usage: samples(n=2, end='se', extra=[])")

def data2run(string):
    """
    data2run :: String -> String
    Return `run` of `data` label
    Example:
    data2run('t2_r1_R1_genes',) -> 't2_r1_genes
    """
    return re.sub(r'_(?:R1|R2)(?=_|$)', '', string)

def data(end='se'):
    """
    data :: Bool -> [(String, Path)]
    Generates iterator of data paths as defined in `sample2list`
    - type='SE' will filter for SE  -> [(key1, Path(R1)) ...]
    - type='PE' will filter for PE  -> [(key11, Path(R1)),(key12, Path(R2)) ...]
    - type='all' will include both SE and PE -> [(key1, Path(R1)),(key11, Path(R1)),(key12, Path(R2)) ...]

    Example:
    data() -> [('t1_r1', PosixPath('../data/test1_R1.fq.gz')),
               ('t2_r1', PosixPath('../data/test2_R1.fq.gz')),
               ('c1_r1', PosixPath('../data/test3_R1.fq.gz'))]

    data(end='pe') ->
              [('t2_r1_R1', PosixPath('../data/SRR23538294_1.fq.gz')),
               ('t2_r1_R2', PosixPath('../data/SRR23538294_2.fq.gz')),
               ('t2_r2_R1', PosixPath('../data/SRR23538293_1.fq.gz')),
               ('t2_r2_R2', PosixPath('../data/SRR23538293_2.fq.gz'))]

    """
    if not end in ['se','pe','all'] : return("usage: data(end='se')")
    if end == 'se':
        return  [(f'{s}_{r}', d[0]) for s,r,d in sample2list if len(d)==1]
    if end == 'pe':
        result = [((f'{s}_{r}_R1',d[0]), (f'{s}_{r}_R2',d[1])) for s,r,d in sample2list if len(d)==2]
        return[(k,r) for k,r in [r for rs in result for r in rs]]
    if end == 'all':
        return data(end='se') + data(end='pe')
    return("usage: data(end='se')")




# def samples2(n=2, se=True, keys=False, extra=[]):
#     '''
#     samples :: Int -> Bool -> Bool -> [String] ->
#     Generates iterators of sample names as defined in 'sample2data'
#     - n=1 returns an iterator of biological replicates, sample(2) -> ['t1', 't2', 'c1']
#     - n=2 (default) returns an iterator of biologial replicates and runs -> ['t1_r1', 't2_r1', 'c1_r1']
#     - se=True (default) filters for SE reads, se=False filters for PE reads
#     - keys=True returns an iterator of tuples of length 'n'. When iterated will require 'n' keys -> [('t1', 'r1'), ('t2', 'r1'), ('c1', 'r1')]
#     - extra=[alist] returns an iterator which is the Cartesian product of the elements of alist and the original iterator
#
#     Example:
#     samples(se=True, keys=True, extra=['gene','genome']) ->
#      [('t1', 'r1', 'gene'), ('t1', 'r1', 'genome'), ('t2', 'r1', 'gene'), ('t2', 'r1', 'genome'), ('c1', 'r1', 'gene'), ('c1', 'r1', 'genome')]
#     '''
#     if n not in [1, 2]: return("usage: samples(n=2, keys=False, extra=False)")
#     if not isinstance(se,bool) or not isinstance(keys,bool): return("usage: samples(n=2, keys=False, extra=False)")
#     valid_extra = isinstance(extra, list) and len(extra) > 0
#     # Generators
#     # [n,keys,extra]
#     results = {
#         1: {(False, False): (s for s in sample2data.keys()),
#             (True, False):  (s for s in sample2data.keys()),
#             (False, True):  (f"{s}_{e}" for s in sample2data.keys() for e in extra),
#             (True, True):   ((s, e) for s in sample2data.keys() for e in extra)},
#
#         2: {(False, False): (f"{s}_{r}" for s, v in sample2data.items() for r in v.keys()),
#             (True, False):  ((s, r) for s, v in sample2data.items() for r in v.keys()),
#             (False, True):  (f"{s}_{r}_{e}" for s, v in sample2data.items() for r in v.keys() for e in extra),
#             (True, True):   ((s, r, e) for s, v in sample2data.items() for r in v.keys() for e in extra)}}
#     gen = results[n].get((keys, valid_extra), "usage: samples(n=2, keys=False, extra=False)")
#
#     if n ==1:
#         return gen
#     # Masks to filter SE and PE Reads
#     if n ==2:
#         # [se][extra]
#         results = {
#             True:  {True:   [len(sample2data[s][r])==1 for s in sample2data for r in sample2data[s]],
#                     False:  [len(sample2data[s][r])==1 for s in sample2data for r in sample2data[s] for e in extra]},
#             False: {True:   [len(sample2data[s][r])==2 for s in sample2data for r in sample2data[s]],
#                     False:  [len(sample2data[s][r])==2 for s in sample2data for r in sample2data[s] for e in extra]}}
#
#         mask = results[se][not extra]
#
#         if se: return compress(gen,mask)
#         if keys: return ((s,r,'R1','R2') for s,r in compress(gen,mask))
#         return ((f'{k}, ''R1','R2') for k in compress(gen,mask))



# def data2(se=True, keys=False):
#     """
#     data :: Bool -> Bool -> [(String, Path)] | [(String, Path, Path)]
#     Generates iterator of data paths as defined in `sample2data`
#     - type='SE' will flatten the list and return a key[s] with a single path to R1 sequence reads
#     - type='PE' will flatten the list and return a key[s] with a two paths to R1,R2  sequence reads

#     Example:
#     data() -> [('t1_r1', PosixPath('../data/test1_R1.fq.gz')),
#                ('t2_r1', PosixPath('../data/test2_R1.fq.gz')),
#                ('c1_r1', PosixPath('../data/test3_R1.fq.gz'))]

#     data(se=False) ->
#               [('t1_r1', PosixPath('../data/test1.1_R1.fq.gz'), PosixPath('../data/test1.1_R2.fq.gz')),
#                ('t1_r2', PosixPath('../data/test_1.2_R1.fq.gz'), PosixPath('../data/test1.2_R2.fq.gz')),
#                ('t1_r3', PosixPath('../data/test1.3_R1.fq.gz'), PosixPath('../data/test1.3_R2.fq.gz'))]
#     """
#     # if type == 'SE': return ((f"{s}_{r}",d) for s, v in sample2data.items() for r, ds in v.items() for d in ds.values())
#     # if type == 'PE': return ((f"{s}_{r}", ds) for s,v in sample2data.items() for r,ds in v.items())
#     # return "usage: data(type), type = {'SE' | 'PE'}"
#     if not isinstance(se,bool): return("usage: data(se=True), for PE  data(se=False)")
#     # [se][keys]
#     results = {
#         True:  {False: ((f'{s}_{r}',sample2data[s][r]['R1']) for s,r in samples2(se=True,keys=True)),
#                 True:  ((s,r,sample2data[s][r]['R1']) for s,r in samples2(se=True,keys=True))},
#         False: {False: ((f'{s}_{r}', sample2data[s][r][r1], sample2data[s][r][r2]) for s,r,r1,r2 in samples2(se=False,keys=True)),
#                 True:  ((s,r, sample2data[s][r][r1], sample2data[s][r][r2]) for s,r,r1,r2 in samples2(se=False,keys=True))}}
#     return results[se][keys]
