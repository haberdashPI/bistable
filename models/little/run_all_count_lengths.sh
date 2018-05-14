N = `cat count_lengths_N.txt`
for i in `seq 1 10 $N`
  sbatch run_count_lengths.sh $((i)) $((i+100)) \
    -r 25 -c 25 --scale_start -1 --scale_top 2 --scale_N 12 \
    -d ~/work/dlittle/bistable_threshold_001/data/ \
    -l ~/work/dlittle/bistable_threshold_001/logs/result_$((i)).log
done

