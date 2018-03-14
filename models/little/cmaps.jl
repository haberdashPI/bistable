import PerceptualColourMaps: cmap
import Colors: RGB
import Colors
using JLD

jldopen(joinpath("packages","AuditoryModel","data","colormaps.jld"),"w") do file
  write(file,"D1",RGB.(cmap("D1")))
  write(file,"reds",RGB.(Colors.colormap("Reds")))
  write(file,"C6",RGB.(cmap("C6")))
end
