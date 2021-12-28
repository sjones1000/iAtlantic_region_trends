#!/bin/bash
#SBATCH --job-name=job_extr_1m_iAtlReg
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=64000
#SBATCH --time=01:00:00
#SBATCH --output=job_extr_1m_iAtlReg.out
#SBATCH --error=job_extr_1m_iAtlReg.err
#SBATCH --partition=cluster

export OMP_NUM_THREADS=14

#
date; set +x

#### latest cdftools
module load matlab_geomar/2020b 

dir_here=/gxfs_work1/geomar/smomw294/iAtlantic/Subprojects/iAtlantic_regions/data/D1p2/VIKING20X_JRA_OMIP/gamma_n/
global_dir=/gxfs_work1/geomar/smomw294/iAtlantic/Subprojects/iAtlantic_regions/data/D1p2/
sim_id=VIKING20X_JRA_OMIP/
in_dir=${global_dir}${sim_id}monthly/
cd $in_dir
in_files=`ls *.nc`
cd $dir_here
out_dir=${global_dir}${sim_id}gamma_n/

#in_files=1_VIKING20X.L46-KFS003_1m_19920101_19921231_grid_T_reg8.nc

for file in ${in_files}; do
	#echo $file
	inpath=\'$in_dir\'
	filen=\'$file\'
	outpath=\'$out_dir\'
	#matlab -nodisplay -nodesktop -r "p1_calc_Yn_worker($in_dir,$file,$out_dir);"
	matlab -nodisplay -nodesktop -nosplash -nojvm -r "try p1_calc_Yn_worker($inpath,$filen,$outpath); catch; end; quit"
done

