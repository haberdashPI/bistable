#!/bin/sh
# S=1
S=101 # start from where N=100 leaves off
# N=100
N=`cat projects/bistable/models/little/count_lengths/count_lengths_N.txt`
# N=2000 # start with just 10 jobs for now, and see how that goes.
K=10

cd projects/bistable
GIT_HASH=`git rev-parse HEAD`
cd

# this just echos the commands, once you verify that it's right, pipe it to sh
for i in `seq $S $K $N`; do
  echo "sbatch projects/bistable/models/little/count_lengths/run_count_lengths.sh $i \
    $((i+K-1)) -r 20 -c 100 --git_hash $GIT_HASH \
    --params projects/bistable/models/little/count_lengths/params_2018-07-16.feather \
    --settings projects/bistable/models/little/count_lengths/settings_2018-07-02.toml \
    -d /scratch/groups/melhila1/dlittle/bistable_threshold_001/data/ \
    -l /scratch/groups/melhila1/dlittle/bistable_threshold_001/logs/result_${i}.log"
done

