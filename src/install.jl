# run this to get started with a newly cloned repository
using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
Pkg.instantiate()

using TOML

config_file = joinpath(@__DIR__,"..","Config.toml")
if isfile(config_file)
  datadir = TOML.parsefile(config_file)["data"]
else
  error("""
Missing a file named `Config.toml` with the following contents:

data = "[folder where bistable data is located]"

This data folder can be found on the lab fileserver under David Little's personal data, under
`data/bistability`. Copy this folder to your own machine and indicate it's
location in `Config.toml`.
""")
end

datadir_link = joinpath(@__DIR__,"..","data")
if !isdir(datadir_link)
  symlink(datadir,datadir_link)
  @info "The folder `data` now links to $datadir"
else
  @info "Directory `data` has already been created."
end
