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
2. Create a file called `Config.toml` in the base directory of the project that
contains the line `data = "[data dir]"` where `[data dir]` is replaced with the location
of the experimental data. This should be the location of the uncompressed archived data.
3. Make sure the julia binary [can be found on your PATH](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started). 
4. Run the install.sh (Mac or Unix) script or install.cmd (Windows). Or you can
run the VSCode build task named "Install".

## Re-running experiment simulations

Installation is all that's necessary to re-run the computational simulations on,
e.g. a cluster. Each experiment, listed under `data`, is run using
the `src/run_all_count_lengths.sh` script.

It is recommended that before running a simulation you should first precompile.

```julia
(bistable> pkg> precompile
```

This will avoid race conditions that can occur in julia's package manager
when running multiple instances of the script across a cluster.

## Interactive use and figure creation

For interactive use and to re-create figures, you should first [install
R](https://cloud.r-project.org/). Then, you'll also want to add the `IJulia`,
`RCall` and `Gadfly`.

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

To reproduce the plot which make use of statistics functions available in R
(Figure 2) you will need to run R and call the following command.
```R
> install.packages('logKDE')
```
