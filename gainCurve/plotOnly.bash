#!/bin/bash

fdr='segbug'
picformat='png'
cp fo.m plotGainCurve.m plotOnly.slurm $fdr
cd $fdr

export picformat
sbatch --export=ALL plotOnly.slurm
