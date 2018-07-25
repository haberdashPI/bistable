#!/bin/sh
#SBATCH
#SBATCH --job-name=bistable
#SBATCH --time=24:0:0
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --requeue
#SBATCH --cpus-per-task=8

module load julia
julia projects/bistable/models/little/count_lengths/run_count_lengths.jl \
	$@ &> work/dlittle/bistable_threshold_002/logs/output_${1}.log

# 40 seconds to run one parameter for a stimulus of with 10 repeats
# that's 40*(50/10) = 200 seconds to run the full length 50 repeat stimulus
# to run 10 repeats for 10 parameters = 5.56 hours
# which I can do for about 100 jobs at once
# so about 5.56 hours x (648/100) = 35 hours
