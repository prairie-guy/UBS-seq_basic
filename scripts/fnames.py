#
# `fnames.py`
#
# C. Bryan Daniels, cdaniels@nandor.net
# 03/17/2024

"""
These functions are generally designed to work with existing `files` and `paths` (in contrast to working
with sample names generated from`config.yaml` . See `configure.py)

As general principles:
    - Use `samples` functions to define file names based upon actual samples as defined in `config.yaml`

    - Use  `fnames` functions when iterating through existing files and paths.

The primary constructors are:
    - `mkpath(path)` and
    - 'fname(path,stem,suffix)`

The primary functions are:
    - `fnames(dir_or_files, suffix="", key='full')`
    - `fnames_string(dir_or_files, delimiter=' ')`
    - `fname_index(path_or_file)`
    - `fname_stem(path_or_file, idx=[0])`
    - `fnames_match(files, pattern, **wildcards)`

Example of a  general workflow in the pipline might be:
  - [pwd() -> workspace]

  - in_path, outpath  = mkpath('path_1'), makepath('path_2'): Create dir withing [workspace]


  - for sample in samples():
      ! <bash command> -i {fnames(in_path,sample,'sam')} -o {fname(out_path, sample, 'bam')}
  or:
  - for sample in fnames(in_path, 'sam'): # By default only the stem is returned
      ! <bash command> -i {fname(in_path,sample,'sam')} -o {fname(out_path, sample, 'bam')}


  - in_path, outpath  = mkpath('path_2'), makepath('path_3'): Create dir withing [workspace]
"""


import os, re
from itertools import product
from subprocess import getoutput as run
from pathlib import Path
from collections.abc import Iterator
from operator import itemgetter as items

def fname(path, stem, suffix, extra=''):
    """
    fname :: str  -> str -> str -> Path
    fname :: Path -> str -> str -> Path

    Return a new Path comprised of a `path`, `stem` and `suffix`
    `extra` appends str to end of stem. Default: extra=''

    Example:
    fname(in_path, 't1_r1', 'bam')          -> Path('some_path/t1_r1.bam)
    fname(in_path, 't1_r1', 'bam', 'genes') -> Path('some_path/t1_r1_genes.bam)
    """

    return Path(path)/f"{stem}.{suffix}" if not extra else Path(path)/f"{stem}_{extra}.{suffix}"

def mkpath(path):
    """
    mkpath :: str -> IO Path

    Creates and returns (not over-writting) a new dir `path`. The current date/time is sent to stdout.

    Example:
    mkpath('dir_name') -> Path('some_path') and new dir `dir_name`

    """
    path = Path(path)
    if not os.path.exists(path): os.makedirs(path)
    date = run('date "+%H:%M:%S_%m-%d-%Y"')
    print(f">>> {{{path}}} {date}")
    return path


def fnames(dir_or_files, suffix="", key='stem'):
    """
    fnames :: str      -> str -> str -> [str]
    fnames :: Path     -> str -> str -> [str]
    fnames :: [str]    -> str -> str -> [str]
    fnames :: Iterator -> str -> str -> [str]

    Returns a list of filenames
    `dir_or_files` can be: a string or Path name of a dir; a list of file names; an Iterator of file names.
    `suffix` is a string by which to filter files, returning only those ending in .`suffix`
    `key` specifies which value of `fname_index(f)` is returned for each file name. Example use of `key`
        fname_index('map_se/t1_r1_genes.bam') ->
            'full':   map_se/t1_r1_genes.bam
            'path':   map_se
            'name':   t1_r1_genes.bam
            'stem':   t1_r1_genes
            'suffix': bam
            'keys':   ['t1', 'r1', 'genes']

    Examples:
    fnames('map_se') ->
        ['map_se/t1_r1_genome.summary', 'map_se/t1_r1_genes.bam.csi', 'map_se/t1_r2_genome.bam', ...]
    fnames('map_se', 'bam') ->
        ['map_se/t1_r2_genome.bam', 'map_se/t1_r2_genes.bam', 'map_se/t1_r1_genome.bam', ...]
    fnames('map_se', 'bam', key='name') ->
        ['t1_r2_genome.bam', 't1_r2_genes.bam', 't1_r1_genome.bam', 'c1_r1_genes.bam', ...]
    fnames('map_se', 'bam', key='stem') ->
        ['t1_r2_genome', 't1_r2_genes', 't1_r1_genome', 'c1_r1_genes', 'c1_r1_genome', 't1_r1_genes']
    fnames('map_se', 'bam', key='keys') ->
        [['t1', 'r2', 'genome'], ['t1', 'r2', 'genes'], ['t1', 'r1', 'genome'], ['c1', 'r1', 'genes'], ...]
    """
    if isinstance(dir_or_files, (str,Path)):
        if not suffix: suffix = '*'
        files = Path(dir_or_files).glob(f'*.{suffix}')
        files = list(map(str,files))

    if isinstance(dir_or_files, list):
        if not suffix: suffix = '*'
        files = list(map(str,dir_or_files))
        files = list(filter(lambda f: re.search(f'\.{suffix}$',f), files)) # check for extension

    if isinstance(dir_or_files, Iterator):
        return fnames(list(dir_or_files), suffix, key) # Remember to update signatue

    return list(map(lambda f: fname_index(f)[key], files))


def fname_index(path_or_file):
    """
    fname_index :: str -> dict
    fname_index :: Path -> dict

    Returns a dict of a file name split by keys=['full', 'path', 'name', 'stem', 'suffix', 'keys'])
    `path_or_file is a string or Path file name`

    Example:
        fname_index('map_se/t1_r1_genes.bam') ->
            'full':   map_se/t1_r1_genes.bam
            'path':   map_se
            'name':   t1_r1_genes.bam
            'stem':   t1_r1_genes
            'suffix': bam
            'keys':   ['t1', 'r1', 'genes']
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


def fnames_string(dir_or_files, delimiter=' '):
    """
    fnames_string :: str      -> str -> str -> [str]
    fnames_string :: Path     -> str -> str -> [str]
    fnames_string :: [str]    -> str -> str -> [str]
    fnames_string :: Iterator -> str -> str -> [str]

    Returns a string concatenated from file names delimited by `delimiter`
    """
    files = fnames(dir_or_files,key='full')
    return delimiter.join(files)


def fnames_match(files, pattern, **wildcards):
    """
    fname_match :: [str] -> str -> [(str, [str])] -> [[str]]
    Returns a list of list of files matching the pattern.

    `files` is a list of file names
    `pattern` is a string similiar to the target `files` with each wildcard indicated by `{ }`
    `pattern` will generally hava a form: {sample}_{run}_{pair}_{ref}.{ext}, though order
    does not matter. Names do not matter provided they correspond to the names used in
    wildcard named lists.
    `**wilcards` are named lists, where the wildcard name is used within the pattern
    For example: sample = ['t1','c1']

    When files are matched to a pattern, chars between two {wildcard} will be ignored.
    This allows grouping of matched files, which are returned as a list

    Example:
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
        combinations = product(*wildcard_values)
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
