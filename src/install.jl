# run this to get started with a newly cloned repository
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using TOML

config_file = joinpath(@__DIR__,"..","Config.toml")
if isfile(config_file)
  datadir = TOML.parsefile(config_file)["data"]
else
  error("""
Missing a file named `Config.toml` with the following contents:

data = "[folder where bistable data is located]"

This folder can be found in the lab fileserver under David Little's personal
data, under `data/bistability`. Copy this file to a locally accessable folder
and point to it in `Config.toml`.
""")
end

symlink(datadir,joinpath(@__DIR__,"..","data"))
