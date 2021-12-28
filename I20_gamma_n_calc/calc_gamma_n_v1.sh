#!/bin/bash
#SBATCH --job-name=job_extr_1m_iAtlReg
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=64000
#SBATCH --time=02:00:00
#SBATCH --output=job_extr_1m_iAtlReg.out
#SBATCH --error=job_extr_1m_iAtlReg.err
#SBATCH --partition=cluster

export OMP_NUM_THREADS=4

#
date; set +x

#### latest cdftools
module load matlab_geomar/2020b 

matlab -nodisplay -nodesktop -nosplash -r "run Calc_Yn_wrapper_script_RUNME; quit"

