#!/usr/bin/env python3

import os, sys, subprocess, resource
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from IPython import get_ipython
from pathlib import Path


# Use iPython command to run shell commands within python
shell = get_ipython().system
get_shell = get_ipython().getoutput

def fname(path, base, sufix):
    'Return a path and suffix complete filename'
    return Path(path)/f"{base}.{sufix}"

def mkpath(path):
    'Return dir path name; creates (not over-writting) a new dir within the pwd'
    path = Path(path)
    if not os.path.exists(path): os.makedirs(path)
    return path
