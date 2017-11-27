# include("units.jl")
include("tempc.jl")
include("stim.jl")

setup_sound(sample_rate=8kHz)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25)
cort = CorticalModel(spect,rates=[2,8,32,64],scales=[0.125,1,4,16])
tempc = TCAnalysis(cort=cort,sparsity=0.8,frame_len=1,prior=8.0,τ=500ms)

# x = playable(sound("../../test.wav"))[0s .. 1s,:left]
# y = cort(spect(x));
# plot_cort(cort,y)

x = aba(60ms,60ms,10,1kHz,6);
# y = cort(spect(x));
# yc = tempc(y);

# rplot(spect,spect(x))
# R"quartz()"
# rplot(tempc,yc)

# better_tempc = TCAnalysis(cort=cort,sparsity=0.1,frame_len=1,prior=8.0,τ=5s)
# byc = better_tempc(spect(x))

# rplot(spect,spect(x))
# R"quartz()"
# rplot(better_tempc,byc)

# TODO: create some graphs of PCA, IPCA, and CCIPCA for cortical
# and spectral representations (also plot spectral and cortical representations
# of the aba paradigm)
# TODO: plotf these representations wrt different f, Δf and Δt

# TODO: setup TC analysis to use online PCA (allowing method to change)
# TODO: plot component representation in spectral and cortical representation

# TODO: implement adaptation and mutual inhibition with a location
# based MI for different layers??
# TODO: determine how to report the object count???

# THOUGHT: we could implement online algorithm by using CCIPCA, with occasional
# refreshes using SVD using the last second or so.

y = spect(x)
a,Y = eigs(y'*y,nev=25)
pca = y*Y

sv, = svds(y,nsv=25)
pca = y*sv[:V]

y = cort(spect(x))
k = OnlinePCA(prod(size(y,2,3,4)),25)
pca = mapslices(y,[2,3,4]) do y_i
  global k = LinAlg.lowrankupdate(k,vec(y_i))
  eigvecs(k)' * vec(y_i)
end

N = 20
y = cort(spect(x))
ii = collect(CartesianRange(size(y,2:4...)))[:];
sv, = svds(y[1:N,ii],nsv=10)
k = OnlinePCA(sv,N,method=:ccipca)
pca = mapslices(y[N:end,:,:,:],[2,3,4]) do y_i
  global k = LinAlg.lowrankupdate(k,vec(y_i))
  eigvecs(k)' * vec(y_i)
end


y = cort(spect(x))
sv, = svds(y[:,ii],nsv=25)
pca = y[:,ii]*sv[:V]

# N = 20
# y = cort(spect(x))
# indices = collect(CartesianRange(size(y,2:4...)))[:]
# sv, = svds(y[1:N,indices],nsv=10)
# k = OnlinePCA(sv,N)
# pca = mapslices(y,[2,3,4]) do y_i
#   global k = push!(k,vec(y_i))
#   eigvecs(k)' * vec(y_i)
# end

# what if we do it all offline?
y = cort(spect(x))
indices = collect(CartesianRange(size(y,2:4...)))[:]
sv, = svds(y[:,indices],nsv=10)



# TODO: figure out what parts of the analysis
# are screwing up the output
#
# leaky covariance (ccipca can fix)
# sparsity (ccipca makes irrelevant)
# cortical representation
#   - note: probably not, since tempc(spect(x)) doesn't work
