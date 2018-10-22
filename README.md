# Bistability Experiment 
WIP experiments to test behavior of the model defined by
[AuditoryBistabilityLE](https://github.com/haberdashPI/AuditoryBistabilityLE).
These files aren't ready for public consumption: there is only minimal effort to
document this code or make it clear to anyone but me, as the approach and
organization may still change dramatically.

# Intallation

These are the steps to initialize this code on a new machine. Mostly for my own reference and
other members of the lab.

1. Install [julia](https://julialang.org/downloads/) version v1.0
2. Open julia in the terminal, from the root directory of this project
3. Run the following commands in the `Pkg` repl (by hitting `]` after opening julia)

```julia
(v1.0) pkg> activate .
(bistability) pkg> add https://github.com/haberdashPI/ShammaModel.jl
(bistability) pkg> add https://github.com/JuliaAudio/SampledSignals.jl#master
(bistability) pkg> add https://github.com/haberdashPI/AuditoryBistabilityLE
(bistability) pkg> instantiate
```

This is all that's necessary to re-run the computational simulations on,
e.g. a cluster. Each experiment, listed under `data`, is run using
the `src/run_all_count_lengths.sh` script.

For interactive use, and to view the results, you'll also want to 
add the `IJulia`, `VegaLite`, and `RCall` pacakges:

```julia
(bistability) pkg> add IJulia VegaLite RCall
```

You can open a Jupyter notebook server and view the notebooks
(under `notebooks` folder) by returning to the julia prompt (by hitting
backspace on an empty line), and calling the `notebook` function
of the `IJulia` package.

```julia
julia> using IJulia; notebook()
```