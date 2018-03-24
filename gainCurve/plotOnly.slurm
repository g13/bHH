#!/bin/bash
#SBATCH -p serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=00:30:00
#SBATCH -o fo_%J.out
#SBATCH --job-name=fO
#SBATCH --mail-user=wd554@nyu.edu
#SBATCH --mail-type=END

# **** Put all #SBATCH directives above this line! ****
# **** Otherwise they will not be in effective! ****
#
# **** Actual commands start here ****
# Load modules here (safety measure)
module purge
module load matlab
ulimit -s unlimited
date

figSave=false
plotTrials=true

matlab -nosplash -r "fo('$picformat','test.cfg',$figSave,$plotTrials);plotGainCurve('$picformat','test.cfg');exit;"
#matlab -nosplash -r "plotGainCurve('$picformat','test.cfg');exit;"

date
hostname
