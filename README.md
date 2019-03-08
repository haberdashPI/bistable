# Bistability Experiment

WIP experiments to test behavior of the model defined by
[AuditoryBistabilityLE](https://github.com/haberdashPI/AuditoryBistabilityLE).
These files aren't ready for public consumption: there is only minimal effort to
document this code or make it clear to anyone but me, as the approach and
organization may still change dramatically.

## Installation

These are the steps to initialize this code on a new machine. Mostly for my own
reference and other members of the lab.

1. Install [Julia](https://julialang.org/downloads/) version v1.0
2. Create a file called `Config.toml` in the base directory of the project and
a line with `data = "[data dir]"` with `[data dir]` replaced with the location
of the experimental data. This can be found the location of the experimental
data. This folder can be found in the lab file server under David Little's
personal data, under `data/bistability`. Copy this file to a locally
accessable folder and point to it in `Config.toml`.
2. Make sure the julia binary [can be found on your PATH](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started). 
2. Run the install.sh (Mac or Unix) script or install.cmd (Windows)

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
