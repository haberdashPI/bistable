using ConjugatePriors
using Distributions
using StatsFuns

@with_kw struct CoherenceTrack
end

const max_sources = 4

struct Weighted{T}
  val::T
  log_weight::Float64
end
sample(x::Weighted{T},x,args...) where T =
  update(sample(x.val,x,args...),x.weight)
update(x::Weighted{T},log_weight::Float64) where T =
  Weighted(x,log_weight+x.log_weight)
log_weight(x::Weighted{T}) = x.log_weight
reweight(x::Weighted{T},log_weight::Float64) = Weighted(x.val,log_weight)
reweight(fn::Function,x) = reweight(x,fn(log_weight(x)))

map_weights(fn::Function,xs) = map(x -> reweight(fn,x),xs)
map_weights!(fn::Function,xs) = map!(x -> reweight(fn,x),xs,xs)

function normalize_weights(x)
  s = logsumexp(log_weight.(x))
  reweight.(w -> w - s,x)
end

resample(N) = xs -> resample(xs,N)
function resample(xs,N)
  weights = rand(Multinomial(N,exp.(log_weight.(xs))))
  filter!(!isinf ∘ log_weight,reweight.(xs,log.(weights)))
end

struct OldParticle{T}
  ϵ::Float64

  z::Vector{T}
  ζ::Vector{Float64}
  v::Float64

  x::Matrix{T}
end
function OldParticle(x::TrackParticle{T})
  OldParticle(p.ϵ,p.z,p.ζ,p.v,p.x)
end

struct TrackParticle{T}
  ϵ::Float64

  z::Vector{Bool}
  ζ::Vector{Float64}
  v::Vector{T}

  x::Matrix{T}
  Σ::Array{Cholesky{T}}

  old::OldParticle{T}
end

function TrackParticle(sources::MVector{SufficientStats},
                       counts::MVector{Float64},
                       active::MVector{Bool},
                       order::MVector{Int},
                       old::TrackParticle{T})
  TrackParticle(sources,counts,active,order,old.n + 1.0,
                old.rate_μ,old.rate_σ,old.scale_μ,old.scale_σ)
end

update(stats::SufficientStats,x::Array{T}) where T = update(stas,vec(x))
update(stats::S,x::Vector{T}) where {T,S} =
  error("No update defined for type $S")
update(stats::MvNormalKnownCovStats{T},x::Vector{T}) where T =
  MvNormalKnownCovStats(stats.invΣ,stats.sx + x,stats.tw + 1.0)

function update(stats::MvNormalStats{T},x::Vector{T}) where T
  s = stats.s + x
  tw = stats.tw + 1.0
  m = s/tw

  z = (s-m)
  s2 = z*z'

  MvNormalStats(s,m,s2,tw)
end

possible_active_sets(n) = []

# TODO: generate list of possible compositions of proposals,
# expressly exclude compositions that leave any components
# that are below some threshold.
function sample(s::TrackState{T},p::TrackPraticle{T},
                proposals::Matrix{T}) where T
  samples = TrackParticle{T}[]

  for r in 1:s.particle_N
    K = length(p.z)

    ϵ = if rand(Bool)
      rand(Normal(0.0,s.σ))
    else
      rand(Normal(p.ϵ,s.ϵ_jump))
    end

    v = if rand(Bool)
      rand(InverseGamma(s.V_α,s.V_θ),K)
    else
      # likelihood-ish update
      abs.(rand(Normal(p.v,s.V_jump),K))
    end
    ζ = rand.(Beta.(p.ζ .* v))
    z = [ζ[k] > rand() for k in 1:K]

    for mapping in possible_mappings(proposals,z)
      x = squeeze(sum(proposals[:,mapping],3),3)
      Σ = similar(p.Σ)
      for i in size(x,2)
        Σ[i] = LinAlg.lowrankupdate(p.Σ[i],x[:,i])
      end

      push!(samples,TrackParticle(ϵ,z,ζ,v,x,Σ,p))
    end
  end

  samples
end

function particle_logpdf(s::TrackState{T},p::TrackParticle{T})
  lp = 0.0

  lp += logsumexp(logpdf(Normal(0.0,s.σ),p.ϵ),
                  logpdf(Normal(p.old.ϵ,s.ϵ_jump),p.ϵ))
  lp += logsumexp([logpdf.(InverseGamma(s.V_α,s.V_θ),p.v);
                   logpdf.(Normal(p.old.v,s.V_jump),p.v)])
  lp += sum(logpdf.(Beta.(p.old.ζ .* p.v),p.ζ))
  lp += sum(p.z ? log(p.ζ[k]) : log(1-p.ζ[k]) for k in 1:k)

  lp
end

struct TrackState{T,M}
  method::M
  particles::Array{Weighted{TrackParticle{T}}}

  # TODO: fill in parameters
end

function logpdf_active(active::Array{Int},counts::AbstractVector{Float64},
                       n::Float64)
  logsumexp(logpdf(i ∈ active ? log(counts/n) : log((1-counts)/n)
                   for i in 1:length(counts)))
end

function logpdf(s::TrackState{T},p::TrackParticle{T},x::Array{T}) where T
  # TODO: implement

  lp = 0.0
  # the observation is a noisy view of a sum of sources
  lp += logpdf(vec(x),MvNormal(squeeze(sum(p.x[:,p.z],2),2),I*p.ϵ))
  lp += logpdf(p.ϵ,Normal(0.0,s.σ))

  # the sources are some variation of their prior values
  lp += logsumexp(
    logpdf.(MvNormal.(p.old.x[:,p.z],PDMats.(p.Σ[p.z])),p.x[:,p.z]))

  # which sources are active is some variation of their prior values
  lp += sum(p.z[k] ? log(p.ζ[k]) : log(1-p.ζ[k]) for k in eachindex(p.z)))
  lp += sum(logpdf.(p.ζ,Beta.(p.old.ζ .* p.v)))
  lp += logpdf(p.v,InverseGamma(s.V_α,s.V_θ))

  lp
end

# basic steps
# generate proposals by each method

# methods:
# draws from the prior
# transformations of old particles
#    max(0,x-y), max(0,y-x)
# draws from the observed data (i.e. no split)
#

# TODO:

# resample proposals, but without removing
# general classes of methods, i.e. maintain
# a range of filters on scales and rates,
# maintain some particles with multiple sources
# maintain some particles from each method???

# for resampling we need to account for the chance of drawing each type of
# particle, each method should assign a probability to each of its own
# proposals, and to any of the other methods proposals
# and we take a weighted sum these probabilities accorording
# to the weight for each method.


# TODO: we can have multiple, parallel particle filters, each with its own
# method

function sample(state::TrackState,y)
  # generate proposals for the components by various methods
  particles = map(weight(state,y),proposals(state,y))
  TrackState(particles |> normalize_weights |> resample(N),state)
end

weight(state::TrackState,y) = p -> weight(state,p,y)
weight(state::TrackState,particle::TrackParticle,y) =
  Weighted(particle,logpdf(state,particle,y) - proposal_logpdf(state,particle))

function sample(states::Array{TrackState{T}},y) where T
  all_particles = Array{Array{Weighted{TrackParticle{T}}}}(length(states))
  for (i,state) in enumerate(states)
    all_particles[i] = map(weight(state,y),proposals(state,y))
    map_weights!(all_particles[i],w -> w - log(length(all_particles[i])))
  end

  log_total = map(all_particles) do particles
    logsumexp(log_weight.(particles) .- log(length(particles)))
  end |> logsumexp

  for particles in all_particles
    map_weights!(particles,w -> w - log_total)
  end

  map(enumerate(states)) do (i,state)
    N = state.particles
    TrackState(all_particles[i] |> resample(N),state)
  end
end

struct FactorizingSampler end

struct NMFObservation
  nmf::NMF.Result
  data
end
components(x::NMFObservation) = x.nmf.W
# TODO: normalize components to ensure sum is accuate
data(x::NMFObservation) = mean(x.data,1)

function proposals(state::TrackState{T,FactorizingSampler} where T)
  [proposal for proposal in sample(state,particle,y.data,components(result))
            for particle in state.particles]
end

function (cohere::CoherenceModel{CoherenceTrack})(x)
  # TODO: apply the particle filter above to each observed step of NMF
  # NOTE: this is a
end
