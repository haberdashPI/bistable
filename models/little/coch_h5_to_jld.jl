push!(LOAD_PATH,"packages")
using AuditoryModel
using JLD
using HDF5

cochba = h5open(joinpath("packages","AuditoryModel","data","cochba.h5")) do file
  read(file,"/real") + read(file,"/imag")*im
end

M = size(cochba,2)
filters = AuditoryModel.CochFilter[]
for ch in 1:M
  p  = floor(Int,real(cochba[1, ch]))
  B  = real(cochba[(0:p)+2, ch])
  A  = imag(cochba[(0:p)+2, ch])

  push!(filters,AuditoryModel.CochFilter(B,A))
end

ch_norm  = imag(cochba[1, M])

data = AuditoryModel.CochFilters(ch_norm,filters)
save(joinpath("packages","AuditoryModel","data","cochba.jld"),"cochba",data)
