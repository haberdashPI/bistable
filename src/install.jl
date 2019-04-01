# run this to get started with a newly cloned repository
using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
Pkg.instantiate()

using TOML

# load the config file
config_file = joinpath(@__DIR__,"..","Config.toml")
if isfile(config_file)
  config = TOML.parsefile(config_file)
  if "data" ∉ keys(config)
    error("Config file missing value for `data`.")
  end
  datadir = config["data"]
else
  error("""
Missing a file named `Config.toml` with the following contents:

data = "[folder where bistable data is located]"

This data folder can be found on the lab fileserver under David Little's
personal data, under `data/bistability`. Copy this folder to your own machine
and indicate it's location in `Config.toml`.

""")
end

# create data link
link = joinpath(@__DIR__,"..","data")
if !isdir(link)
  symlink(datadir,link)
  @info "The folder `$link` now links to $datadir"
else
  @info "Directory `$link` has already been created."
end

# generate the workspace file
str = """
{
	"folders": [
		{
			"path": "."
		},
		{
			"path": "$(joinpath(Pkg.devdir(),"AuditoryBistabilityLE"))"
		},
		{
			"path": "$(joinpath(Pkg.devdir(),"ShammaModel.jl"))"
		},
		{
			"path": "$(joinpath(Pkg.devdir(),"PlotAxes.jl"))"
		},
		{
			"path": "$(joinpath(Pkg.devdir(),"MetaArrays.jl"))"
		}
	],
	"settings": {},
	"extensions": {
		"recommendations": [
			"bungcip.better-toml",
			"ikuyadeu.r",
			"julialang.language-julia",
			"colinfang.markdown-julia",
			"davidanson.vscode-markdownlint",
			"haberdashPI.terminal-polyglot"
		]
	}
}
"""
open("bistable.code-workspace",write=true) do f
  write(f,str)
end
