# Bistability Experiment

Experiments to test behavior of the model defined by
[AuditoryBistabilityLE](https://github.com/haberdashPI/AuditoryBistabilityLE).

Part of the in-submission paper:

Little, D.F. Snyder, J.S., Elhilali, M., Ensemble modeling of auditory
streaming reveals potential sources of bistability across the perceptual
hierarchy.

### Abstract 

Perceptual bistability---the spontaneous fluctuation of perception between
two interpretations of a stimulus---occurs when observing a large variety of
ambiguous stimulus configurations. This phenomenon has the potential to serve
as a tool for, among other things, understanding how function varies across
individuals due to the large individual differences that manifest during
perceptual bistability. Yet it remains difficult to interpret the functional
processes at work, without knowing where bistability arises during
perception. In this study we explore the hypothesis that bistability
originates from multiple sources distributed across the perceptual hierarchy.
We develop a hierarchical model of auditory processing comprised of three
distinct levels: a Peripheral, tonotopic analysis, a Central analysis
computing features found more centrally in the auditory system, and an Object
analysis, where sounds are segmented into different streams. We model
bistable perception within this system by injecting adaptation, inhibition
and noise into one or all of the three levels of the hierarchy. We evaluate a
large ensemble of variations of this hierarchical model, where each model has
a different configuration of adaptation, inhibition and noise. This approach
avoids the assumption that a single configuration must be invoked to explain
the data. Each model is evaluated based on its ability to replicate two
hallmarks of bistability during auditory streaming: the selectivity of
bistability to specific stimulus configurations, and the characteristic
log-normal pattern of perceptual switches. Consistent with a distributed
origin, a broad range of model parameters across this hierarchy lead to a
plausible form of perceptual bistability. The ensemble also appears to
predict that greater individual variation in adaptation and inhibition occurs
in later stages of perceptual processing.

## Installation

These are the steps to initialize this code on a new machine.

1. Install [Julia](https://julialang.org/downloads/) version v1.1
2. Create a file called `Config.toml` in the base directory of the project that
contains the line `data = "[data dir]"` where `[data dir]` is replaced with the location
of the experimental data. This should be the location of the uncompressed archived data.
You can use the same data as reported in our in-submission paper by downloading
this [archive](https://osf.io/9shxv/) of our data
3. Make sure the julia binary [can be found on your PATH](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started). 
4. Run the install.sh (Mac or Unix) script or install.cmd (Windows). Or you can
run the VSCode build task named "Install".

### Developer tools

The project is well set up to be easily edited in [Visual Studio
Code](https://code.visualstudio.com/). When you open
`bistable.code-workspace` inside VSCode (a file created by `install.sh`), the
project comes with recommendations for all of the extensions you will need to
read and write code for the project. You can also search and edit all of the
supporting packages written for this project within VSCode as well (e.g.
AuditoryBistabilityLE).

## Code Organization

There are a few principles for code organization:

Interactive code---material that involved examining results and frequent
changes to code based on this output---is stored in Jupyter notebooks,
organized by the figure numbers from the submitted paper.

Implementation of the specific simulations is organized under `src`. It contains
all of the material that is specific to the particular simulations and analyses
run in the paper.

General purpose code---computational models or general-use data
structures---are organized into packages. These will be installed in
`$JULIA_PKG_DEVDIR` (default location is `~/.julia/dev/`) when you call
`install.sh` or `install.cmd`. The packages are:

* `AuditoryBistabilityLE` this is all of the code specific to the model
implemented in the paper. It is described by the `Model design` subsection of
the paper.
* `ShammaModel.jl`. Implements the auditory spectrogram; a high-level port of
  the basic functionality in the NSLtoolbox for matlab.
* `PlotAxes.jl`. General utility for creating quick plots in Julia.
* `MetaArrays.jl`. General utility for adding meta-data to arrays in Julia.

For more details please refer to the README.md files of the individual packages.

The latter three packages will eventually be published to Julia public
repositories, as they have more general-purpose uses (particularly the latter
two), but were implemented in the process of running these computational
simulations.

## Re-running computational simulations

After installing, each experiment, listed under `data`, can be run using the
`src/run_all_count_lengths.sh` script. This uses
[SLURM](https://slurm.schedmd.com/documentation.html) to run and manage jobs
on a computer cluster.

It is recommended that before running a simulation you should first precompile.

```julia
(bistable> pkg> precompile
```

This will avoid race conditions that can occur in julia's package manager
when running multiple instances of the script across a cluster.

## Interactive use and figure creation

After re-running simulations, or using the previously run simulations, you
can interact with the results. For this interactive use, and to re-create
figures, you should first [install R](https://cloud.r-project.org/). Then,
you'll also want to add the `IJulia`, `RCall`, `Fontconfig`, `Cairo` and `Gadfly`.

```julia
(bistability) pkg> add IJulia RCall Fontconfig Cairo Gadfly
```
You can open a Jupyter notebook server and view the notebooks
(under `notebooks` folder) by returning to the Julia prompt (by hitting
backspace on an empty line), and calling the `notebook` function
of the `IJulia` package.

```julia
julia> using IJulia; notebook()
```

Alternatively, you can open the notebook server directly on the command line,
using jupyter.

```bash
$ jupyter notebook
```

Or you can use jupter lab, for a more modern notebook experience.

```bash
$ jupyter lab
```

### Figure 2
To reproduce the plot which makes use of statistics functions available in R
(Figure 2) you will need to run R and call the following command.
```R
> install.packages('logKDE')
```

### Figure reproducibility 

Note that the raw figures will look somewhat different than the final
published form: all of the figures were loaded into Adobe Illustrator and
edited for clarity. (The data themselves were not modified.)
