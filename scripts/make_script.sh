#!/usr/bin/env bash

cd $(dirname "$(readlink -f "$0")")
cd ../workspace/
cp ubs_basic.py ubs_basic.older.py
jupyter nbconvert --to script ubs_basic.ipynb
