#!/bin/sh
#SBATCH
#SBATCH --job-name=d_bistable
#SBATCH --time=0:10:0
#SBATCH --partition=debug
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2

module load julia
JULIA_NUM_THREADS=2 julia projects/bistable/models/little/run_count_lengths.jl \
	$@ &> work/dlittle/bistable_threshold_001/logs/output_$((1)).log
