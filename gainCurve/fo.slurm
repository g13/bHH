#!/bin/bash
#SBATCH -p serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH -o fo_%J.out
#SBATCH --mail-user=wd554@nyu.edu
#SBATCH --mail-type=END

# **** Put all #SBATCH directives above this line! ****
# **** Otherwise they will not be in effective! ****
#
# **** Actual commands start here ****
# Load modules here (safety measure)
module purge
module load gcc
module load matlab
module load boost
date

./fo && matlab -nosplash -r "fo('$picformat','test.cfg');exit;"
date
hostname
