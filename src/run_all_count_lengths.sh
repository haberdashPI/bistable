#!/bin/sh
result_dir=$1
label=$2
repeat=${3:-20}
stim_count=${4:-100}
proj_dir="projects/bistable/src"

S=1
# S=101 # start from where N=100 leaves off
# N=10
N=`cat ${result_dir}/${label}_N.txt`
# N=2000 # start with just 10 jobs for now, and see how that goes.
K=200

cd ${proj_dir}
GIT_HASH=`git rev-parse HEAD`
cd

# this just echos the commands, once you verify that it's right, pipe it to sh
for i in `seq $S $K $N`; do
  echo "sbatch ${proj_dir}/run_count_lengths.sh $i \
    $((i+K-1)) -r ${repeat} -c ${stim_count} --git_hash $GIT_HASH \
    --params ${result_dir}/params.jld2 \
    --settings ${proj_dir}/settings.toml \
    -d ${result_dir}/data/ \
    -l ${result_dir}/logs/result_${i}.log"
done

