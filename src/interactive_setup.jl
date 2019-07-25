include("base_setup.jl")
using RCall
R"options(rcalljl_options=list(width=800,height=400))"
# import Cairo, Fontconfig
using Revise # (we have to call this in IJulia... seems like a bug)
revise()

includet(joinpath(srcdir,"count_lengths.jl"))
includet(joinpath(srcdir,"measures.jl"))
includet(joinpath(srcdir,"plotting.jl"))
