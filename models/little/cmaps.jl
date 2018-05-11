import PerceptualColourMaps: cmap
import Colors: RGB
import Colors
using JLD2

jldopen(joinpath("packages","AuditoryModel","data","colormaps.jld2"),"w") do file
  file["D1"] = RGB.(cmap("D1"))
  file["reds"] = RGB.(Colors.colormap("Reds"))
  file["C6"] = RGB.(cmap("C6"))
end
