#!/bin/bash
#SBATCH -p serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=17
#SBATCH --mem=60G
#SBATCH --time=4:00:00
##SBATCH -a 1-1
##SBATCH -o ReR._%A_%a.out
#SBATCH -o sL_%J.out
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
