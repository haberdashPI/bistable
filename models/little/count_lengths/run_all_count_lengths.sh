#!/bin/sh
S=1
# S=2201 # start from where N=2000 leaves off
N=`cat projects/bistable/models/little/count_lengths_N.txt`
# N=2000 # start with just 10 jobs for now, and see how that goes.
K=200
for i in `seq $S $K $N`; do
  echo "sbatch projects/bistable/models/little/run_count_lengths.sh $i \
    $((i+K-1)) -r 10 -c 25 --scale_start -1 --scale_stop 2 --scale_N 15 \
    -d /scratch/groups/melhila1/dlittle/bistable_threshold_001/data/ \
    -l /scratch/groups/melhila1/dlittle/bistable_threshold_001/logs/result_${i}.log"
done

