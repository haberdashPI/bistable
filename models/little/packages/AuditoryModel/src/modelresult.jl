using AxisArrays

abstract type ModelResult{T,N} <: AbstractArray{T,N} end
# ModelResult children should implement 2 methods
# data(x) = an AxisArray of the data
# params(x) = an immutable struct of parameters for the given model
# optionally you can implement resultname(x) to make printing prettier

Base.size(x::ModelResult) = size(data(x))
@inline @Base.propagate_inbounds Base.getindex(x::ModelResult,i...) =
  modelwrap(x,getindex(data(x),i...))
@inline @Base.propagate_inbounds Base.getindex(x::ModelResult,i::Int...) =
  getindex(data(x),i...)
@inline @Base.propagate_inbounds Base.setindex!(x::ModelResult,v,i...) =
  setindex!(data(x),v,i...)
@inline @Base.propagate_inbounds Base.setindex!(x::ModelResult{T},v::T,
                                                i::Int...) where T =
  setindex!(data(x),v,i...)

modelwrap(x::ModelResult{T},val::T) where T = val
modelwrap(x::A,val::AxisArray{T,N}) where {T,N,A <: ModelResult{T,N}} =
  A(val,params(x))
modelwrap(x::A,val::AxisArray{T}) where {T,A <: ModelResult{T}} = val
modelwrap(x::ModelResult,val::ModelResult) = val
modelwrap(x::ModelResult,val::T) where T =
  error("unexpected result type $T from getindex on a `ModelResult`.")

Base.IndexStyle(x::ModelResult) = IndexStyle(data(x))
function Base.similar(x::ModelResult,::Type{S},
                      dims::NTuple{N,Int}) where {S,N}
  if length(dims) == ndims(x)
    typeof(x)(similar(data(x),S,dims),x.params)
  else
    similar(data(x),S,dims)
  end
end

function showparams(io,x)
  for field in fieldnames(x)
    write(io,"$field = ")
    show(io,getfield(x,field))
    write(io,"\n")
  end
end

resultname(x::ModelResult) = string(typeof(x))

function Base.show(io::IO,x::ModelResult{T,N}) where {T,N}
  write(io,resultname(x))
  write(io,":\n")
  showparams(io,params(x))
  write(io,"-------\n")
  write(io,"$(N)d $T data with axes: \n")
  for ax in axes(data(x))
    write(io,"$(AxisArrays.axisname(ax)): $(round(ustrip(ax.val[1]),2))"*
          " - $(round(ustrip(ax.val[end]),2)) $(unit(ax.val[1]))\n")
  end
end
Base.show(io::IO,::MIME"text/plain",x::ModelResult) = show(io,x)

AxisArrays.axisdim(x::ModelResult,ax) = axisdim(data(x),ax)
AxisArrays.axes(x::ModelResult,i...) = axes(data(x),i...)
