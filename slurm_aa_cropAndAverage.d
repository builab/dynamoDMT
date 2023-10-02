#!/bin/bash
#SBATCH --ntasks=1                
#SBATCH --partition=cpu
#SBATCH --job-name=crop # Job name
#SBATCH --error=aa_cropnaverage.err
#SBATCH --output=aa_cropnaverage.out
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=15GB

module load matlab
module load cuda-8.0.61-gcc-4.9.4-rv3d2jh


matlab -nodisplay < dynamoDMT/aa_cropAndAverage.m
