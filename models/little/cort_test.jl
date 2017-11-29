# include("units.jl")
include("tempc.jl")
include("stim.jl")

setup_sound(sample_rate=8kHz)

spect = AuditorySpectrogram("/Users/davidlittle/Data/cochba.h5",len=25)
cort = CorticalModel(spect,rates=[2,8,32,64],scales=[0.125,1,4,16])

# TODO: setup TC analysis to use online or offline PCA (allowing method to
# change in constructor)

# TODO: plotf these representations wrt different f, Δf and Δt

# TODO: implement adaptation and mutual inhibition with a location
# based MI for different layers??
# TODO: determine how to report the object count???

# TODO: plot component representation in spectral and cortical representation?

#=========================================
# iterate through various pca algorithms
dir = "../../plots/run_2017_11_28"
x = aba(60ms,60ms,10,1kHz,6);
params = [(:pca,0),(:ccipca,0),(:ipca,0),(:ccipca,20),(:ipca,20)]
for (represent,str) in [(cort,"cort")] # (spect,"spect")
  y = represent(x)

  rplot(represent,y)
  R"ggsave(paste($dir,'/',$str,'_aba.pdf',sep=''),width=8,height=6)"
  for (method,init_len) in params
    tempc = TCAnalysis(represent,10,method=method,init_len=init_len)
    tc = tempc(y)

    rplot(tempc,tc)
R"""
    ggsave(paste($dir,'/',$str,'_',$(string(method)),'_init',
           $(string(init_len)),'.pdf',sep=''),width=8,height=6)
"""

    rplot(tempc,abs.(tc))
R"""
    ggsave(paste($dir,'/',$str,'_',$(string(method)),'_abs_init',
           $(string(init_len)),'.pdf',sep=''),width=8,height=6)
"""

  end
end
=========================================#

# cleanup of iteratation through various algorithms
x = aba(60ms,60ms,10,1kHz,6);
dir = "../../plots/run_2017_11_28_PCA_algorithms"
params = [(:pca,0),(:ipca,0),(:ccipca,20)]
for (represent,str,width,height) in [(cort,"cort",10,6),(spect,"spect",8,6)]
  y = represent(x)

  rplot(represent,y)
  R"ggsave(paste($dir,'/',$str,'_aba.pdf',sep=''),width=$width,height=$height)"
  for (method,init_len) in params
    tempc = TCAnalysis(represent,10,method=method,init_len=init_len)
    tc = tempc(y)

    rplot(tempc,abs.(tc))
R"""
    ggsave(paste($dir,'/',$str,'_',$(string(method)),'_abs_init',
           $(string(init_len)),'.pdf',sep=''),width=8,height=6)
"""

  end
end

# vary the stimulus parameters
dir = "../../plots/run_2017_11_28_Stimulus_parameters"
function plot_stages(x,dir,prefix)
  for (represent,str,width,height) in [(cort,"cort",10,6),(spect,"spect",8,6)]
    y = represent(x)

    rplot(represent,y)
R"""
    ggsave(paste($dir,'/',$prefix,$str,'.pdf',sep=''),
           width=$width,height=$height)
"""

    tempc = TCAnalysis(represent,10,method=:ipca,init_len=0)
    tc = tempc(y)

    rplot(tempc,abs.(tc))
R"""
    ggsave(paste($dir,'/',$prefix,$str,'_pca.pdf',sep=''),width=8,height=6)
"""

  end
end

for delta in [1,6,12]
  x = aba(60ms,60ms,10,1kHz,delta)
  plot_stages(x,dir,"delta$(delta)_")
end

for freq in [500Hz,1000Hz,2000Hz]
  x = aba(60ms,60ms,10,freq,6)

  plot_stages(x,dir,"freq$(ustrip(freq))Hz_")
end

for len in [30ms,42ms,60ms,85ms,120ms,240ms]
  x = aba(len,len,10,1kHz,6)

  plot_stages(x,dir,"length$(ustrip(len))ms_")
end
