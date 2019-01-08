# Bistability Experiment

WIP experiments to test behavior of the model defined by
[AuditoryBistabilityLE](https://github.com/haberdashPI/AuditoryBistabilityLE).
These files aren't ready for public consumption: there is only minimal effort to
document this code or make it clear to anyone but me, as the approach and
organization may still change dramatically.

## Installation

These are the steps to initialize this code on a new machine. Mostly for my own
reference and other members of the lab.

1.  Install [Julia](https://julialang.org/downloads/) version v1.0
2.  Open Julia in the terminal, from the root directory of this project
3.  Run the following commands in the `Pkg` REPL (by hitting `]` after opening
    Julia)

```julia
(v1.0) pkg> activate .
(bistable) pkg> add https://github.com/wildart/TOML.jl#v0.4.0
(bistable) pkg> add https://github.com/haberdashPI/MetaArrays.jl
(bistable) pkg> add https://github.com/JuliaAudio/SampledSignals.jl#master
(bistable) pkg> add https://github.com/haberdashPI/ShammaModel.jl
(bistable) pkg> add https://github.com/haberdashPI/AuditoryBistabilityLE
(bistable) pkg> instantiate
```

This is all that's necessary to re-run the computational simulations on,
e.g. a cluster. Each experiment, listed under `data`, is run using
the `src/run_all_count_lengths.sh` script.

It is recommended that before running a simulation you should first precompile.

```julia
(bistable> pkg> precompile
```

This will avoid race conditions that can occur in julia's package manager
when running multiple instances of the script across a cluster.

## Interactive use and figure creation

For interactive use, and to recreate the figures, you'll also want to
add the `IJulia`, `RCall` and `Gadfly`. You should [install
R](https://cloud.r-project.org/) first.

```julia
(bistability) pkg> add IJulia Gadfly RCall
```
You can open a Jupyter notebook server and view the notebooks
(under `notebooks` folder) by returning to the Julia prompt (by hitting
backspace on an empty line), and calling the `notebook` function
of the `IJulia` package.

```julia
julia> using IJulia; notebook()
```

To reproduce the plots which make use of statistics functions available in R
(Figure 3) you will need to run R and call the following command.
```R
> install.packages('logKDE')
```
## Experimental Notes

TODO: describe format of jld2 files (they are hdf5 files)
