using NMF
using Parameters

@with_kw struct CoherenceNMF
  normalize::Bool = true
  maxiter::Int = 2000
  tol::Float64 = 1e-4
end
CoherenceMethod(::Type{Val{:nmf}},params) = CoherenceNMF(;params...)

Base.eltype(::CoherenceNMF,::Type) = Float64
Base.eltype(::CoherenceNMF,::AbstractArray) = Float64

function with_method(foreach_window,method::CoherenceNMF,K)
  convergence_count,total_count = 0,0
  # Winit,Hinit = fill(0.0,(0,0)),fill(0.0,(0,0))
  # method = NMF.ALSPGrad{Float64}(tol=tcparams.method.tol,
  #                                maxiter=tcparams.method.maxiter)

  foreach_window() do window
    x_t = reshape(window,prod(size(window,1,2)),:)
    k = min(size(x_t,1),K)

    # if isempty(Winit)
    #   Winit,Hinit = NMF.nndsvd(abs.(window),k)
    # end

    # solution = NMF.solve!(method,abs.(window),Winit,Hinit)
    total_count += 1

    solution = nnmf(abs.(x_t) .+ method.tol/2,k,init=:nndsvd,
                    tol=method.tol,
                    maxiter=method.maxiter)

    convergence_count += !solution.converged
    solution.converged || warn("NMF failed to converge.",once=true,
                               key=object_id(method))

    y = reshape(solution.H,:,size(window)[3:end]...)
    permutedims(y,[2:ndims(window)-1;1])
  end

  if convergence_count > 0
    percent = 100round(convergence_count / total_count)
    if percent > 0
      info("$percent% of frames failed to fully converge to a solution.")
    else
      info("<0.01% of frames failed to fully converge to a solution.")
    end
  end
end
