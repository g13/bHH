#!/bin/bash
#SBATCH -p serial
#SBATCH --nodes=1
#SBATCH --ntasks=8
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --mem=40G
#SBATCH --mem=16G
##SBATCH --time=00:05:00
#SBATCH --time=04:00:00
#SBATCH -a 1-1
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
theme='test_k4'
model='HH'
picformat='png'
draw=true
#draw=false
ppp=true
#ppp=false
#load=true
load=false
#name='FS_somato_Rat'
name='RS_exc_Rat'
npool=$SLURM_NTASKS
echo $npool
echo $SLURM_TASKS_PER_NODDE
group=(1)
#group=(1 2 3 4 5 6 7 8 9 10 11 12 13 14)
#i=1
#group=(7)
echo ${#group[@]}
echo $SLURM_ARRAY_TASK_MAX
date
if [ ${#group[@]} -eq $SLURM_ARRAY_TASK_MAX ]
then
    i=$SLURM_ARRAY_TASK_ID
    pick=${group[$i-1]} 
    echo $i
    echo $pick
    
    ############################################### CHECK YOUR WALL TIME ##################################################################
    #matlab -nosplash -nodesktop -nodisplay -r "noAdapV('$theme','$name',$pick,'$model','$picformat',$draw,$ppp,$load,$npool);exit;"
    matlab -nosplash -nodesktop -nodisplay -r "getLib;exit;"
    #matlab -nosplash -nodesktop -nodisplay -r "single_check_Dev('$name',$pick,'$model','$picformat',$draw,$ppp,$load);exit;"
fi
date
hostname
