#!/bin/sh
S=1
label=$1
# S=101 # start from where N=100 leaves off
# N=10
N=`cat projects/bistable/src/${label}_count_lengths_N.txt`
# N=2000 # start with just 10 jobs for now, and see how that goes.
K=10

cd projects/bistable
GIT_HASH=`git rev-parse HEAD`
cd

proj_dir="projects/bistable/src"
result_dir="/scratch/groups/melhila1/dlittle"

# this just echos the commands, once you verify that it's right, pipe it to sh
for i in `seq $S $K $N`; do
  echo "sbatch ${proj_dir}/run_count_lengths.sh $i \
    $((i+K-1)) -r 20 -c 100 --git_hash $GIT_HASH \
    --params ${proj_dir}/${label}_params_2018-09-07.feather \
    --settings ${proj_dir}/settings_2018-09-07.toml \
    -d ${result_dir}/bistable_threshold_${label}/data/ \
    -l ${result_dir}/bistable_threshold_${label}/logs/result_${i}.log \
    ${result_dir}/bistable_threshold_${label}/logs/output_${i}.log"
done

