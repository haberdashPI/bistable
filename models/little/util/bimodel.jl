
include(joinpath(@__DIR__,"biscales.jl"))
using AuditoryCoherence

function build_priors(scales,early,N)
  p = reshape(mean(early,4),size(early,1),:)
  base_prior = AuditoryCoherence.IsoMultiNormalStats(p,N);

  map(scales) do c
    cur_prior = deepcopy(base_prior)
    cur_prior.S *= c;
    cur_prior
  end
end

function bistable_model(cs,params,args;progressbar=false,
                        intermediate_results=false)
  csa = bistable_scales(cs,params,args,progressbar=progressbar,
                        intermediate_results=intermediate_results)
  rates = 2.^linspace(args["rate_start"],args["rate_stop"],
                      args["rate_N"])
  crs = cortical(csa[:,:,args["min_Hz"] .. resoultion["max_Hz"]];
                 rates=[-rates;rates])
  C = cohere(
    crs,
    ncomponents=args["nmf_K"],
    window=args["nmf_window"]*ms,
    method=:nmf,
    delta=args["nmf_delta"]*ms,
    maxitor=args["nmf_itr"],
    tol=1e-3,
    progressbar=progressbar
  )

  freq_N = args["track_prior_freq_N"]
  freq_bias = args["track_prior_freq_bias"]

  Ct,lp,tc,ps = track(
    C,
    method=:multi_prior,
    max_sources = args["track_max_sources"],
    tcs = linspace(args["track_tc_start"],args["track_tc_stop"],
                   args["track_N"]),
    thresh=args["track_zero_thresh"],

    source_priors = build_priors([1.0],C[0s .. args["track_early_C"]*s],
                                 args["track_prior_strength"])
    freq_prior = AuditoryCoherence.BinomialCond(
      :old => AuditoryCoherence.Beta(freq_N*(1-freq_bias),freq_N),
      :new => AuditoryCoherence.Beta(freq_N*freq_bias,freq_N)
    ),
    progressbar=progressbar
  )

  lengths = count_streams(Ct,lp)

  if intermediate_results
    lengths,csa,C,Ct,lp,tc,ps
  else
    lengths
  end
end





