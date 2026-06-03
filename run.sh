#!/bin/bash

# initialize Conda for this script's shell
source ~/miniconda3/etc/profile.d/conda.sh

# Activate env
conda activate snakemake

# Run snakemake 
snakemake --cores all --sdm conda --rerun-incomplete
