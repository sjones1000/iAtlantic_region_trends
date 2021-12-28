#!/bin/bash
#SBATCH --job-name=gn_reg6
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=64000
#SBATCH --time=08:00:00
#SBATCH --output=gn_reg6.out
#SBATCH --error=gn_reg6.err
#SBATCH --partition=cluster

export OMP_NUM_THREADS=4

#
date; set +x

#### latest cdftools
module load matlab_geomar/2020b 

matlab -nodisplay -nodesktop -nosplash -r "run Calc_Yn_wrapper_script_RUNME_reg6; quit"

