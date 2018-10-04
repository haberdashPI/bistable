#!/bin/sh
dir=$1
label=$2
proj_dir="projects/bistable/src"
result_dir="/scratch/groups/melhila1/dlittle/${dir}"

S=1
# S=101 # start from where N=100 leaves off
# N=10
N=`cat ${result_dir}/${label}_N.txt`
# N=2000 # start with just 10 jobs for now, and see how that goes.
K=15

cd ${proj_dir}
GIT_HASH=`git rev-parse HEAD`
cd

# this just echos the commands, once you verify that it's right, pipe it to sh
for i in `seq $S $K $N`; do
  echo "sbatch ${proj_dir}/run_count_lengths.sh $i \
    $((i+K-1)) -r 20 -c 100 --git_hash $GIT_HASH \
    --params ${result_dir}/${label}_params.feather \
    --settings ${proj_dir}/settings.toml \
    -d ${result_dir}/data/ \
    -l ${result_dir}/logs/result_${i}.log"
done

