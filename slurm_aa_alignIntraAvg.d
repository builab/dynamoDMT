#!/bin/bash
#SBATCH --ntasks=1                
#SBATCH --partition=titan
#SBATCH --job-name=alignIntraAvg # Job name
#SBATCH --error=aa_alignIntraAvg.err
#SBATCH --output=aa_alignIntraAvg.out
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:6
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=15GB

module load matlab
module load cuda-8.0.61-gcc-4.9.4-rv3d2jh


matlab -nodisplay < dynamoDMT/aa_alignIntraAvg.m

