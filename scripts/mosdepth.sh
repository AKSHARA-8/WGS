#!/usr/bin/bash

bamfile=$1
output_prefix=$2

source activate new_base

mosdepth -m -n -t 3 --fast-mode --by 500 ${output_prefix} ${bamfile}

conda deactivate
