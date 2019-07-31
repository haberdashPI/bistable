# run this to get started with a newly cloned repository
using Pkg
root = joinpath(@__DIR__,"..")
Pkg.activate(root)
Pkg.Registry.add(RegistrySpec(url = "https://github.com/haberdashPI/LabRegistry.jl"))

# we use a base.Manifest file, so that we can specify a local directory for
# the AuditoryBistabilityLE package.  This allows the user to edit this package
# as desired.
cp(joinpath(root,"base.Manifest.toml"),joinpath(root,"Manifest.toml"))
Pkg.instantiate()
Pkg.develop("AuditoryBistabilityLE")

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

This data folder can be downloaded from https://osf.io/9shxv/ and stored
almost anywhere you like; it should not be stored inside this project
directory.

""") end

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
