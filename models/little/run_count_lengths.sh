#!/bin/sh
#SBATCH
#SBATCH --job-name=bistable
#SBATCH --time=6:0:0
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --requeue
#SBATCH --cpus-per-task=8

module load julia
JULIA_NUM_THREADS=8 julia projects/bistable/models/little/run_count_lengths.jl \
	$@ &> work/dlittle/bistable_threshold_001/logs/output_${1}.log

# 37 seconds for 4 parameters, in 2 threads with 2 repeats and 10 ab's
# to run a single parameter set for 100 repeats with 25 abs' with
# 16 threads should take 37/2/2/10*25*100/16/60 = 2.40885 minutes
# to run 100 of these parameters in a single job should take
# 37/2/2/10*25*100/16/60*100/60 = 4.015 hours.
