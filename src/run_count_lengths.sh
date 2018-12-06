#!/bin/sh
#SBATCH
#SBATCH --job-name=bistable
#SBATCH --time=0:20:0
#SBATCH --partition=debug
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --requeue
#SBATCH --cpus-per-task=2

# module load julia
julia  -e 'using Pkg; Pkg.activate("projects/bistable")' \
       -e 'include("projects/bistable/src/run_count_lengths.jl")' \
       -O3 --banner=no $@

