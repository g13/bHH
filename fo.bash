#!/bin/bash
dir=$1
if [ "$fdr" == "" ]; then
    echo using default folder "test"
    dir='test'
fi
cd gainCurve
if [ -d "$dir" ]; then
    rm -r $dir/*
else
    mkdir $dir
fi
cp fo $dir
cp read_cfg.m $dir
cp fo.m $dir
cp fo.slurm $dir
cp ../input.cfg $dir/test.cfg

cd $dir
sbatch fo.slurm
