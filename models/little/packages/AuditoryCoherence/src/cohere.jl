using Sounds
using MacroTools
using Match
using Parameters
using Parameters
using AxisArrays
using Unitful

export cohere, component, mask, ncomponents, components, component_means,
  windowing, map_windowing

abstract type CoherenceMethod end
struct CParams{M,P} <: AuditoryModel.Params
  cort::P
  ncomponents::Int
  skipframes::Int
  window::typeof(1.0s)
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
ncomponents(x::AxisArray) = size(x,axisdim(x,Axis{:component}))
components(x::CParams) = 1:ncomponents(x)
components(x::AuditoryModel.Result) = components(AxisArray(x))
components(x::AxisArray) = axisvalues(axes(x,Axis{:component}))[1]

component(x::Union{AxisArray,Coherence},n) = x[Axis{:component}(n)]

AuditoryModel.hastimes(x::Coherence) = HasTimes()

function component_means(C)
  mdims = filter(x -> x != axisdim(C,Axis{:component}),1:ndims(C))
  vec(mean(C,mdims))
end

AuditoryModel.frame_length(params::CParams,x) =
  max(1,floor(Int,params.delta / Δt(x)))

function CParams(x;ncomponents=1,window=1s,
                  method=:nmf,delta=10ms,skipframes=0,
                  normalize_phase=true,method_kwds...)
  method = CoherenceMethod(Val{method},method_kwds)

  CParams(AuditoryModel.Params(x), ncomponents, skipframes,
          convert(typeof(1.0s),window), convert(typeof(1.0s),delta), method)
end

windowing(x,params::CParams) =
  windowing(x,length=params.window,step=params.delta)

windowing(x,dim=timedim(x);kwds...) = windowing(hastimes(x),x,dim;kwds...)
map_windowing(fn,x,dim=timedim(x);kwds...) =
  map_windowing(fn,hastimes(x),x,dim;kwds...)
function map_windowing(fn,::HasTimes,x,dim;step=nothing,kwds...)
  windows = windowing(x,dim;step=step,kwds...)
  xs = map(windows) do ixs
    fn(x[Axis{:time}(ixs)])
  end
  AxisArray(xs,axes(windows,Axis{:time}))
end
map_windowing(fn,::HasNoTimes,x,dim;kwds...) =
  map(ixs -> fn(x[Axis{:time}(ixs)]),windowing(HasNoTimes(),x,dim;kwds...))

function windowing(::HasNoTimes,x::AbstractArray,dim;
                   length=nothing,step=nothing,minlength=length)
  (max(1,t-length+1):t for t in indices(x,dim)[minlength:step:end])
end

function windowing(::HasTimes,data::AbstractArray,dim;
                   length=nothing,step=nothing,minlength=length)
  helper(x::Number) = x
  helper(x::Quantity) = max(1,floor(Int,x / Δt(data)))
  length_,step_,minlength_ = helper.((length,step,minlength))

  win = windowing(HasNoTimes(),data,dim,
                  length=length_,step=step_,minlength=minlength_)
  AxisArray(collect(win),
            Axis{:time}((minlength_:step_:size(data,dim))*Δt(data)))
end

windowlen(params::CParams,x) = round(Int,params.window/Δt(x))
function nunits(params::CParams,x)
  mapreduce(*,axes(x)) do ax
    isa(ax,Axis{:time}) || isa(ax,Axis{:rate}) ? 1 : length(ax)
  end
end

cohere(x::AuditoryModel.Result;progressbar=true,params...) =
  cohere(x,CParams(x;params...),progressbar)

function cohere_progress(progressbar,x,params)
  if progressbar
    windows = windowing(x,params)
    Progress(length(windows),desc="Temporal Coherence Analysis: ")
  end
end

function cohere(x::AbstractArray,params::CParams,progressbar=true,
                progress = cohere_progress(progressbar,x,params))
  @assert axisdim(x,Axis{:time}) == 1
  @assert axisdim(x,Axis{:rate}) == 2
  @assert ndims(x) == 4 # AxisArray can't handle skipped indices so we assum
                        # the right dimensionality

  # if we already have components, just wrap up the values with
  # parameters (since we've already computed components)
  if :component in axisnames(x)
    return Coherence(x,params)
  end

  windows = windowing(x,params)

  K = ncomponents(params)
  C_data = zeros(eltype(params.method,x),length(windows),size(x)[3:end]...,K)
  C = AxisArray(C_data,
                axes(windows,Axis{:time}),
                axes(x)[3:end]...,
                Axis{:component}(1:K))

  with_method(params.method,K) do extract
    for (i,w_inds) in enumerate(windows)
      skipped = w_inds[1:1+params.skipframes:end]
      # axis array can't handle skipped indices, so we assume
      # the right dimensionality
      components = extract(x.val.data[skipped,:,:,:])
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

  windows = enumerate(windowing(cr,AuditoryModel.Params(C)))
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

  cortical(y,AuditoryModel.Params(cr))
end
