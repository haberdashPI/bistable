# run this to get started with a newly cloned repository
using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
Pkg.instantiate()

using TOML

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

This data folder can be found on the lab fileserver under David Little's personal data, under
`data/bistability`. Copy this folder to your own machine and indicate it's
location in `Config.toml`.
""")
end

function makelink(dir,link)
  if isnothing(link)
    return
  end

  if !isdir(link)
    symlink(dir,link)
    @info "The folder `$link` now links to $dir"
  else
    @info "Directory `$link` has already been created."
  end
end
makelink(joinpath(@__DIR__,"..","data"),datadir)
makelink(joinpath(@__DIR__,"..","juliadev"),Pkg.devdir())
