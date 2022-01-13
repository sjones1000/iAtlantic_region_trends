#!/bin/bash
#SBATCH --job-name=V20_timeseries
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=64000
#SBATCH --time=2:00:00
#SBATCH --output=V20_timeseries.out
#SBATCH --error=V20_timeseries.err
#SBATCH --partition=cluster

export OMP_NUM_THREADS=4

#
date; set +x

#### latest cdftools
module load matlab_geomar/2020b 

matlab -nodisplay -nodesktop -nosplash -r "run Timeseries_plot_wrapper_script_RUNME.m; quit"

