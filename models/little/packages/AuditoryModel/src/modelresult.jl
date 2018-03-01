using AxisArrays

abstract type Params end
abstract type ModelResult{T,N} <: AbstractArray{T,N} end
# ModelResult children should implement 2 methods
# AxisArray(x) = an AxisArray of the data
# Params(x) = a Params object
# optionally you can implement resultname(x) to make printing prettier

Base.size(x::ModelResult) = size(AxisArray(x))
@inline @Base.propagate_inbounds Base.getindex(x::ModelResult,i...) =
  modelwrap(x,getindex(AxisArray(x),i...))
@inline @Base.propagate_inbounds Base.getindex(x::ModelResult,i::Int...) =
  getindex(AxisArray(x),i...)
@inline @Base.propagate_inbounds Base.setindex!(x::ModelResult,v,i...) =
  setindex!(AxisArray(x),v,i...)
@inline @Base.propagate_inbounds Base.setindex!(x::ModelResult{T},v::T,
                                                i::Int...) where T =
  setindex!(AxisArray(x),v,i...)

modelwrap(x::ModelResult{T},val::T) where T = val
modelwrap(x::A,val::AxisArray{T,N}) where {T,N,A <: ModelResult{T,N}} =
  A(val,Params(x))
modelwrap(x::A,val::AxisArray{T}) where {T,A <: ModelResult{T}} = val
modelwrap(x::ModelResult,val::ModelResult) = val
modelwrap(x::ModelResult,val::T) where T =
  error("unexpected result type $T from getindex on a `ModelResult`.")

Base.IndexStyle(x::ModelResult) = IndexStyle(AxisArray(x))

# TODO: add similar version consistent with
# similar from AxisArrays which takes axes objects instead of dims
function Base.similar(x::ModelResult,::Type{S},
                      dims::NTuple{N,Int}) where {S,N}
  if size(x) == dims
    data = similar(AxisArray(x),S)
    # TODO: I need to figure out
    similar_helper(x,data,x.params)
  else
    similar(AxisArray(x),S,dims)
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
  showparams(io,Params(x))
  write(io,"-------\n")
  write(io,"$(N)d $T data with axes: \n")
  for ax in axes(AxisArray(x))
    write(io,"$(AxisArrays.axisname(ax)): $(round(ustrip(ax.val[1]),2))"*
          " - $(round(ustrip(ax.val[end]),2)) $(unit(ax.val[1]))\n")
  end
end
Base.show(io::IO,::MIME"text/plain",x::ModelResult) = show(io,x)

AxisArrays.axisdim(x::ModelResult,ax) = axisdim(AxisArray(x),ax)
AxisArrays.axes(x::ModelResult,i...) = axes(AxisArray(x),i...)
