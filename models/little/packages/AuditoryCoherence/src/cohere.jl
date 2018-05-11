using Sounds
using MacroTools
using Match
using Parameters
using Parameters
using AxisArrays

export cohere, component, mask, ncomponents, components, component_means

abstract type CoherenceMethod end
struct CParams{M,P} <: AuditoryModel.Params
  cort::P
  ncomponents::Int
  window::typeof(1.0s)
  minwindow::typeof(1.0s)
  delta::typeof(1.0s)
  method::M
end

struct Coherence{M,T,N} <: AuditoryModel.Result{T,N}
  val::AxisArray{T,N}
  params::CParams{M}
end
AuditoryModel.Params(x::Coherence) = x.params
AxisArrays.AxisArray(x::Coherence) = x.val
AuditoryModel.resultname(x::Coherence) = "Coherence Components"
AuditoryModel.similar_helper(::Coherence,val,params) = Coherence(val,params)
AuditoryModel.Δt(as::CParams) = as.delta

struct CoherenceComponent{M,T,N} <: AuditoryModel.Result{T,N}
  val::AxisArray{T,N}
  params::CParams{M}
end
AuditoryModel.Params(x::CoherenceComponent) = x.params
AxisArrays.AxisArray(x::CoherenceComponent) = x.val
AuditoryModel.resultname(x::CoherenceComponent) = "Single Coherence Component"
AuditoryModel.similar_helper(::CoherenceComponent,val,params) =
  CoherenceComponent(val,params)

function AuditoryModel.modelwrap(x::A,newval::AxisArray{T}) where
  {T,M,A <: Coherence{M,T}}

  if setdiff(axisnames(x),axisnames(newval)) == [:component]
    CoherenceComponent(newval,AuditoryModel.Params(x))
  elseif axisnames(x) == axisnames(newval)
    A(newval,AuditoryModel.Params(x))
  else
    newval
  end
end

function Coherence(x::Coherence,p::CParams)
  @assert x.params == p "Coherence parameters do not match"
  x
end

function Coherence(x::AbstractArray,p::CParams)
  @assert nfreqs(p.cort) == nfreqs(x) "Frequency channels do not match"
  Coherence(x,p)
end

ncomponents(x::CParams) = x.ncomponents
ncomponents(x::AuditoryModel.Result) = length(components(x))
components(x::CParams) = 1:ncomponents(x)
components(x::AuditoryModel.Result) = components(AxisArray(x))
components(x::AxisArray) = axisvalues(axes(x,Axis{:component}))[1]

component(x::Coherence,n) = x[Axis{:component}(n)]

function component_means(C)
  mdims = filter(x -> x != axisdim(C,Axis{:component}),1:ndims(C))
  vec(mean(C,mdims))
end

AuditoryModel.frame_length(params::CParams,x) =
  max(1,floor(Int,params.delta / Δt(x)))

function CParams(x;ncomponents=1,window=1s,minwindow=window,
                  method=:nmf,delta=10ms,
                  normalize_phase=true,method_kwds...)
  method = CoherenceMethod(Val{method},method_kwds)

  CParams(AuditoryModel.Params(x),ncomponents,
          convert(typeof(1.0s),window),
          convert(typeof(1.0s),minwindow),
          convert(typeof(1.0s),delta),
          method)
end

windowing(x,dim,params::CParams) =
  windowing(x,dim,length=windowlen(params,x),step=frame_length(params,x),
            minlength=min_windowlen(params,x))
windowing(x,dim;length=nothing,step=nothing,minlength=nothing) =
  (max(1,t-length+1):t for t in indices(x,dim)[minlength:step:end])

windowlen(params::CParams,x) = round(Int,params.window/Δt(x))
min_windowlen(params::CParams,x) = round(Int,params.minwindow / Δt(x))
function nunits(params::CParams,x)
  mapreduce(*,axes(x)) do ax
    isa(ax,Axis{:time}) || isa(ax,Axis{:rate}) ? 1 : length(ax)
  end
end

cohere(x::AuditoryModel.Result;params...) = cohere(x,CParams(x;params...))

function cohere(x::AbstractArray{T},params::CParams) where T
  @assert axisdim(x,Axis{:time}) == 1
  @assert axisdim(x,Axis{:rate}) == 2

  # if we already have components, just wrap up the values with
  # parameters (since we've already computed components)
  if :component in axisnames(x)
    return Coherence(x,params)
  end

  windows = windowing(x,1,params)

  K = ncomponents(params)
  C_data = zeros(eltype(params.method,x),length(windows),size(x)[3:end]...,K)
  C = AxisArray(C_data,
                Axis{:time}(times(x)[map(last,windows)]),
                axes(x)[3:end]...,
                Axis{:component}(1:K))

  with_method(params.method,K) do extract
    progress = Progress(length(windows),desc="Temporal Coherence Analysis: ")
    for (i,w_inds) in enumerate(windows)
      components = extract(x[Axis{:time}(w_inds)])
      C[i,indices(components)...] = components

      next!(progress)
    end
  end

  Coherence(C,params)
end

function mask(cr::AbstractArray,C::Coherence)
  error("Please select one component (see documentation for `component`).")
end

function mask(cr::AbstractArray{T},C::CoherenceComponent) where T
  @assert axisdim(cr,Axis{:time}) == 1
  @assert axisdim(cr,Axis{:rate}) == 2
  @assert size(cr)[3:end] == size(C)[2:end] "Dimension mismatch"

  windows = enumerate(windowing(cr,1,AuditoryModel.Params(C)))
  y = zeros(AxisArray(cr))
  norm = similar(y,real(T))
  norm .= zero(real(T))
  @showprogress "Masking: " for (i,w_inds) in windows
    c = C[Axis{:time}(i)]
    y[Axis{:time}(w_inds)] .+= reshape(c,1,1,size(c)...)
    norm[Axis{:time}(w_inds)] += 1
  end
  y ./= norm
  y ./= maximum(abs,y)
  y .= sqrt.(abs.(cr) .* y) .* exp.(angle.(cr).*im)

  @show nfreqs(y)
  cortical(y,AuditoryModel.Params(cr))
end
