include(joinpath(@__DIR__,"count_lengths.jl"))

# summraize the number of parameters run for each *.clbin file specified
counts = zeros(27000)
dir = "work/dlittle/bistable_threshold_001/data"
for file in readdir(dir)
  if ismatch(r"clbin$",file)
    for row in loadrows(joinpath(dir,file))
      if row.param_index < length(counts)
        counts[row.param_index] += 1
      end
    end
  end
end


