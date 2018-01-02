#!/bin/bash
#SBATCH -p serial
#SBATCH --nodes=1
#SBATCH --ntasks=24
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --mem=40G
#SBATCH --mem=16G
##SBATCH --time=00:05:00
#SBATCH --time=04:00:00
##SBATCH -a 1-1
##SBATCH -o ReR._%A_%a.out
#SBATCH -o EIF.%J.out
#SBATCH --mail-user=wd554@nyu.edu
#SBATCH --mail-type=END

# **** Put all #SBATCH directives above this line! ****
# **** Otherwise they will not be in effective! ****
#
# **** Actual commands start here ****
# Load modules here (safety measure)
module purge
module load matlab
date
    matlab -nosplash -nodesktop -nodisplay -r "getLib;exit;"
date
hostname
