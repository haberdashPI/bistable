#!/bin/sh

result_dir=$1
label=$2
K=${3:-10}
repeat=${4:-20}
default_N=`cat ${result_dir}/${label}_N.txt`
N=${5:-$default_N}

proj_dir="projects/bistable/src"
S=1

# S=101 # start from where N=100 leaves off
# N=10
# N=2000 # start with just 10 jobs for now, and see how that goes.

cd ${proj_dir}
GIT_HASH=`git rev-parse HEAD`
cd

# this just echos the commands, once you verify that it's right, pipe it to sh
for i in `seq $S $K $N`; do
  echo "sbatch ${proj_dir}/run_count_lengths.sh $i \
    $((i+K-1)) -r ${repeat} --git_hash $GIT_HASH \
    --params ${result_dir}/params.jld2 \
    --settings ${proj_dir}/settings.toml \
    -d ${result_dir}/data/ \
    -l ${result_dir}/logs/result_${i}.log"
done

