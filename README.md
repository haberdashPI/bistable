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

For interactive use, and to view the results, you'll also want to
add the `IJulia`, `VegaLite`, and `RCall` packages:

```julia
(bistability) pkg> add IJulia VegaLite RCall
```

You can open a Jupyter notebook server and view the notebooks
(under `notebooks` folder) by returning to the Julia prompt (by hitting
backspace on an empty line), and calling the `notebook` function
of the `IJulia` package.

```julia
julia> using IJulia; notebook()
```
