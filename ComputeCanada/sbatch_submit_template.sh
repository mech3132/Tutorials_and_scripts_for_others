#!/bin/bash

#SBATCH --account=def-USER
#SBATCH --time=35:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=EMAIL@gmail.com
#SBATCH --mail-type=ALL


cd ~/scratch # Change into scratch directory
mkdir -p outputs # Make outputs directory in scratch
### ALL INPUTS SHOULD BE IN SCRATCH FOLDER ####


module load singularity

singularity exec -B $PWD:/home -B ~/scratch:/inputs -B ~/scratch/outputs:/outputs \
 ~/images/qiime2-2021.11.sif \
  qiime fragment-insertion sepp \
  --i-representative-sequences /repset_16S.qza \
  --i-reference-database /sepp-refs-silva-128.qza \
  --o-tree /outputs/insertion-tree_16S.qza \
  --o-placements /outputs/insertion-placements_16S.qza
 
 
