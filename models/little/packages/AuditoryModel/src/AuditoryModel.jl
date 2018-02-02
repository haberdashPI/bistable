module AuditoryModel
using DataFrames

# flag to determine weather this package compiles with or without support for
# MATLAB implementation of functions.
const USING_MATLAB = false

export AuditorySpectrogram, freqs, times, scales, rates, freq_ticks, rplot,
    frame_length, Δt, delta_t, Δf, delta_f, CorticalModel, plot_scales,
    plot_scales2, Seconds, Hertz, TimeDim, FreqDim

include("units.jl")
include("plots.jl")
include("audio_spect.jl")
include("cortical.jl")

end
