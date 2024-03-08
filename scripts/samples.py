#!/usr/bin/env python3

import os, sys
from pathlib import *
from subprocess import getoutput
run = getoutput
from snakemake import load_configfile
from snakemake.io import expand
from collections import defaultdict

home_path      = Path.cwd()/'..'
config = load_configfile(home_path/"config.yaml")
SAMPLES = config['samples']

def fname(path, base, sufix):
    'Return a path and suffix complete filename'
    return Path(path)/f"{base}.{sufix}"

def mkpath(path):
    'Return dir path name; creates (not over-writting) a new dir within the pwd. Also prints date/time executed'
    path = Path(path)
    if not os.path.exists(path): os.makedirs(path)
    date = run('date +%H:%M:%S_%m-%d-%Y')
    print(f">>> {{{path}}} {date[0]}")
    return path

read_ids = ["R1", "R2"]
pairend_run_ids = []
sample2data = defaultdict(dict)
group2sample = defaultdict(list)
for s, v in config["samples"].items():
    if v.get("treated_", True):
        group2sample[v["group"]].append(s)
    for i, v2 in enumerate(v["data"], 1):
        r = f"r{i}"
        sample2data[s][r] = {k3: home_path/v3 for k3, v3 in v2.items()}
        if len(v2) == 2:
            pairend_run_ids.append(s + "_" + r)

[f"{s}_{r}" for s, v in sample2data.items() for r in v.keys()]

[(f"{s}_{r}" ,d) for s, v in sample2data.items() for r, ds in v.items() for d in ds.values()]
[(s,r,ds) for s,v in sample2data.items() for r,ds in v.items()]

def samples(x='samples', keys=False, extra=None):
    "make generators"
    # n = 1
    if x == 'samples' and keys: return ((s,r) for s,v in sample2data.items() for r in v.keys())
    if x == 'samples': return (f"{s}_{r}" for s,v in sample2data.items() for r in v.keys())

    #if x == 'keys': return (s for s in sample2data.keys())
 # n = 2
    if x == 'values': return (v for v in sample2data.values())
    if x == 'items': return ((s,v) for s,v in sample2data.items())

    if x == 'data': return ((s,r,ds) for s,v in sample2data.items() for r,ds in v.items())


def samples(keys=False, extra=False):
    valid_extra    = isinstance(extra, list) and len(extra) > 0
    results = {
        (False, False): (f"{s}_{r}" for s,v in sample2data.items() for r in v.keys()),
        (True, False): ((s,r) for s,v in sample2data.items() for r in v.keys()),
        (False, True): (f"{s}_{r}_{e}" for s,v in sample2data.items() for r in v.keys() for e in extra) ,
        (True, True): ((s,r,e) for s,v in sample2data.items() for r in v.keys() for e in extra) }
    return results.get((keys, valid_extra), "Invalid combination")


def samples(n=2, keys=False, extra=False):
    if n not in [1, 2]: return "Invalid value for n"
    valid_extra = isinstance(extra, list) and len(extra) > 0
    results = {
        1: {(False, False): (s for s in sample2data.keys()),
            (True, False):  (s for s in sample2data.keys()),
            (False, True):  (f"{s}_{e}" for s in sample2data.keys() for e in extra),
            (True, True):   ((s, e) for s in sample2data.keys() for e in extra)},

        2: {(False, False): (f"{s}_{r}" for s, v in sample2data.items() for r in v.keys()),
            (True, False):  ((s, r) for s, v in sample2data.items() for r in v.keys()),
            (False, True):  (f"{s}_{r}_{e}" for s, v in sample2data.items() for r in v.keys() for e in extra),
            (True, True):   ((s, r, e) for s, v in sample2data.items() for r in v.keys() for e in extra)}}

    return results[n].get((keys, valid_extra), "Usage: samples(n=2, keys=False, extra=False)")


# {'t1':
#  {'description': 'HeLa polyA+ RNA treated with ultrafast BS, replicate 3',
#   'group': 'treated',
#   'other': 'SRR23538290',
#   'data': [{'R1': 'data/test1_R1.fq.gz'},
#            {'R1': 'data/test1.1_R1'},
#            {'R1': 'data/test1.2_R1', 'R2': 'data/test1.2_R2'}]},
#  't2':
#  {'description': 'HeLa polyA+ RNA treated with ultrafast BS, replicate 2',
#   'group': 'treated',
#   'other': 'SRR23538291',
#   'data': [{'R1': 'data/test2_R1.fq.gz'}]},
# 'c1':
#  {'description': 'HeLa polyA+ RNA treated with ultrafast BS, replicate 1',
#   'group': 'control',
#   'other': 'SRR23538292',
#   'data': [{'R1': 'data/test3_R1.fq.gz'}]}}

# samples(keys) ->
#
# ['t1','t2','c1']
# list(sample2data.keys())
# samples(runs) ->
# ((s,r) for s,v in sample2data.items() for r in v.keys())
# ['t1_r1', 't1_r2', 't1_r3', 't2_r1', 'c1_r1']

#samples(data, type='SE') ->
# [(f"{s}_{r}",d) for s,v in sample2data.items() for r,ds in v.items() for d in ds.values()]
# [('t1_r1', PosixPath('/home/cdaniels/uofc_data/ubs_seq/UBS-seq_basic/scripts/../data/test1_R1.fq.gz')),
#  ('t1_r2', PosixPath('/home/cdaniels/uofc_data/ubs_seq/UBS-seq_basic/scripts/../data/test1.1_R1')),
#  ('t1_r3', PosixPath('/home/cdaniels/uofc_data/ubs_seq/UBS-seq_basic/scripts/../data/test1.2_R1')),
#  ('t1_r3', PosixPath('/home/cdaniels/uofc_data/ubs_seq/UBS-seq_basic/scripts/../data/test1.2_R2')),
#  ('t2_r1', PosixPath('/home/cdaniels/uofc_data/ubs_seq/UBS-seq_basic/scripts/../data/test2_R1.fq.gz')),
#  ('c1_r1', PosixPath('/home/cdaniels/uofc_data/ubs_seq/UBS-seq_basic/scripts/../data/test3_R1.fq.gz'))]

# samples(data,type='PE')
# ((s,r,ds) for s,v in sample2data.items() for r,ds in v.items())
# [('t1_r1', {'R1': PosixPath('/home/cdaniels/uofc_data/ubs_seq/UBS-seq_basic/scripts/../data/test1_R1.fq.gz')}),
#  ('t1_r2', {'R1': PosixPath('/home/cdaniels/uofc_data/ubs_seq/UBS-seq_basic/scripts/../data/test1.1_R1')}),
#  ('t1_r3', {'R1': PosixPath('/home/cdaniels/uofc_data/ubs_seq/UBS-seq_basic/scripts/../data/test1.2_R1'),
#             'R2': PosixPath('/home/cdaniels/uofc_data/ubs_seq/UBS-seq_basic/scripts/../data/test1.2_R2')}),
#  ('t2_r1', {'R1': PosixPath('/home/cdaniels/uofc_data/ubs_seq/UBS-seq_basic/scripts/../data/test2_R1.fq.gz')}),
#  ('c1_r1', {'R1': PosixPath('/home/cdaniels/uofc_data/ubs_seq/UBS-seq_basic/scripts/../data/test3_R1.fq.gz')})]
