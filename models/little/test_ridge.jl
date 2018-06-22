
push!(LOAD_PATH,"packages")
using AuditoryModel
using AuditoryCoherence
include("util/stim.jl")

x        = ab(120ms,120ms,1,40,500Hz,6) |> normpower |> amplify(-10dB)
sp       = audiospect(x)
cs       = cortical(sp;scales=cycoct.*round.(2.0.^linspace(-2,2.5,9),1))

crs = cortical(cs[:,:,400Hz .. 800Hz];rates=[(-2.0.^(1:5))Hz; (2.0.^(1:5))Hz])
C = cohere(crs,ncomponents=3,window=150ms,method=:nmf,skipframes=2,
           delta=100ms,maxiter=100,tol=1e-3)

prior = ridgenorm(C,1,scale=0.25,freq=0.25)
prior.μ .= 0
alert()

AuditoryCoherence.logpdf(prior,zero(vec(C[2,:,:,1])))
AuditoryCoherence.logpdf(prior,vec(C[1,:,:,1]))

prior2 = isonorm(C,1)
prior.μ .= 0

AuditoryCoherence.logpdf(prior2,zero(vec(C[1,:,:,1])))
AuditoryCoherence.logpdf(prior2,vec(C[1,:,:,1]))

# TODO: okay, we've tested out a few things for logpdf_mvt
# now try looking at the results of tracking again

# The key test will be examining logpdf's with low and high scales
