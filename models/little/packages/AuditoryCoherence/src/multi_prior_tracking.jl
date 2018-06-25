
@with_kw struct MultiPriorTracking <: Tracking
  tcs::Array{typeof(1.0s)} = [100ms, 250ms, 500ms, 1s, 2.5s, 5s]
  thresh::Float64 = 1e-2
  source_priors
  freq_prior
  unmodeled_prior::Float64 = 1e-2
  max_sources::Int = 5
end
Tracking(::Val{:multi_prior};params...) = MultiPriorTracking(;params...)

function expand_params(params::MultiPriorTracking)
  [PriorTracking(tc,params.thresh,prior,params.freq_prior,
                 params.unmodeled_prior,params.max_sources)
   for tc in params.tcs for prior in params.source_priors]
end

function track(C::Coherence,params::MultiPriorTracking,progressbar=true,
               progress=track_progress(progressbar,ntimes(C),"multi-prior"))
  @show params.max_sources
  time = Axis{:time}
  component = Axis{:component}
  # hopefully one day I can remove these assertions
  # but not right now
  @assert axisdim(C,time) == 1
  @assert axisdim(C,component) == 4

  tracks = TrackedSources.(expand_params(params))
  C_out_tr = [similar(C,axes(C)[1:end-1]...,component(1:params.max_sources))
              for i in 1:length(tracks)]
  lp_out = [AxisArray(fill(0.0,ntimes(C)),axes(C,1))
            for i in 1:length(tracks)]
  for i in 1:length(tracks); C_out_tr[i] .= 0; end
  decays = [1.0 - 1.0 / max(1.0,track.params.tc / Δt(C)) for track in tracks]

  @showprogress "Tracking Sources..." for t in eachindex(times(C))
    for (track_index,track) in enumerate(tracks)
      decay = decays[track_index]
      MAPgrouping = maximumby(grouping -> logpdf(track,C[time(t)],grouping),
                              possible_groupings(params.max_sources,
                                                 ncomponents(C)))
      lp_out[track_index][time(t)] = logpdf(track,C[time(t)],MAPgrouping)

      # arrange the sources
      for (kk,i) in iterable(MAPgrouping)
        for k in kk
          #  C_out[time(t),component(i)] .+= C[time(t),component(k)]
          # I hope in the future I can use the above but not now
          # TODO: submit an issue with AxisArrays
          C_out_tr[track_index][t,:,:,i] .+= C[t,:,:,k]
        end
      end

      # update the source models
      mult!(track,decay)
      update!(track,C_out_tr[track_index][time(t)],MAPgrouping)
    end
  end

  for i in eachindex(C_out_tr)
    order = sortperm(component_means(C_out_tr[i]),rev=true)
    C_out_tr[i]  .= C_out_tr[i][component(order)]
  end

  tcs = [track.params.tc for track in tracks]
  priors = [mean(std(track.params.source_prior)) for track in tracks]
  C_out_tr, lp_out, tcs, priors
end
