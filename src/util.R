library(JuliaCall)
julia = julia_setup()
julia_command("using Pkg; Pkg.activate(\"..\")")
julia_command("using AuditoryBistabilityLE")

percept_lengths = function(x,minlen,framerate){
  julia_assign("over",x)
  julia_assign("minlen",floor(minlen / framerate))
  result=julia_eval("AuditoryBistabilityLE.mergelengths(AuditoryBistabilityLE.findlengths(over)...,minlen)")
  data.frame(length=result[[1]]*framerate,stimulus=as.numeric(result[[2]]))
}

interpolate_times = function(x,times,to){
  julia_call("interpolate_times",x,times,to)
}

W = function(x){
  if (sum(!is.na(x)) > 5000){
    shapiro.test(sample(x,5000))[[1]]
  }else if (sum(!is.na(x)) < 3){
    0
  }else if (all(is.na(x)) || all(x == first(x[!is.na(x)]),na.rm=T)){
    0
  }else{
    shapiro.test(x)[[1]]
  }
}

set_bound = function(xs){
  if (length(xs) == 0){
    logical(0)
  } else if (length(xs) == 1){
    c(T)
  }else if (length(xs) == 2){
    c(T,T)
  }else {
    c(T,rep(F,length(xs)-2),T)
  }
}

clamp = function(x,lower,upper){
  pmin(upper,pmax(lower,x))
}

clean_ratio = function(stimulus,length,is_bound){
  clean_length = sum(length[!is_bound])
  full_length = sum(length)

  if (full_length - clean_length > 0.2*full_length){
    stim = stimulus
    len = length
    use_full=T
  }else{
    stim = stimulus[!is_bound]
    len = length[!is_bound]
    use_full=F
  }

  find_ratio = function(stim,len,use_full){
    m1 = mean(log10(len[stim == 0]))
    m2 = mean(log10(len[stim == 1]))
    clamp(m1 - m2,-1,1)
  }

  if (any(stim == 0) && any(stim == 1)){
    find_ratio(stim,len)
  }else{
    stim = stimulus
    len = length
    use_full=T
    if (any(stim == 0) && any(stim == 1)){
      find_ratio(stim,len)
    }else if (any(stimulus == 0)){
      1
    }else{
      -1
    }
  }
}

