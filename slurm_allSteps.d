#!/bin/bash
#SBATCH --ntasks=1                
#SBATCH --partition=titan
#SBATCH --job-name=intraAln # Job name
#SBATCH --error=aa_intraAln_all.err
#SBATCH --output=aa_intraAln_all.out
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:6
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=15GB

module load matlab
module load cuda-8.0.61-gcc-4.9.4-rv3d2jh

# If you are very confident about your data. Not recommended to run everything all at once

matlab -nodisplay < dynamoDMT/aa_cropAndAverage.m
matlab -nodisplay < dynamoDMT/aa_intraAln.m
matlab -nodisplay < dynamoDMT/aa_alignIntraAvg.m
#matlab -nodisplay < dynamoDMT/aa_generateAxonemeAvg.m
matlab -nodisplay < dynamoDMT/aa_alignAllParticles.m
matlab -nodisplay < dynamoDMT/aa_filamentRepick.m
matlab -nodisplay < dynamoDMT/aa_intraAlnRepick.m
matlab -nodisplay < dynamoDMT/aa_alignRepickAvg.m
matlab -nodisplay < dynamoDMT/aa_alignRepickParticles.m





