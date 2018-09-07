using Logging
using Dates

struct DatedLogger <: AbstractLogger
  stream::IOStream
  level::LogLevel
end
DatedLogger(stream::IOStream) = DatedLogger(stream,Logging.Info)
Logging.shouldlog(logger::DatedLogger, level, _module, group, id) = true
Logging.min_enabled_level(logger::DatedLogger) = logger.level

function Logging.handle_message(logger::DatedLogger, level, message,
                                _module, group, id, filepath, line, keys...)

  # copied and modified from
  # https://github.com/JuliaLang/julia/blob/v1.0.0/base/logging.jl

  buf = IOBuffer()
  iob = IOContext(buf, logger.stream)
  levelstr = level == Logging.Warn ? "Warning" : string(level)
  msglines = split(chomp(string(message)), '\n')
  time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

  filepath = if startswith(something(filepath,""),@__DIR__)
    filepath[length(@__DIR__)+2:end]
  elseif occursin(r"AuditoryBistabilityLE/src/",filepath)
    part = match(r"AuditoryBistabilityLE/src/(.+)$",filepath)
    "[bi-model]/"*part[1]
  elseif occursin(r"ShammaModel/src/",filepath)
    part = match(r"ShammaModel/src/(.+)$",filepath)
    "[aud-mode]/"*part[1]
  else
    filepath
  end

  # use a short, single line message if possible
  offset = 60
  if length(msglines) == 1 && length(filepath) < 30 &&
    length(msglines[1]) <= 80

    headerbuf = IOBuffer()
    print(headerbuf, "[ ", levelstr, "[", time, "] "," ",
          something(filepath,"nothing"), ":", something(line,"xx")," ")

    header = String(take!(headerbuf))
    print(iob,header)
    for r in 1:(offset - length(header))
      print(iob," ")
    end
    print(iob,"|")
    println(iob, msglines[1])
  # for longer messages, fall back to the standard long format.
  else
    println(iob, "┌ ", levelstr, "[", time, "]: ", msglines[1])
    for i in 2:length(msglines)
      println(iob, "│ ", msglines[i])
    end
    for (key, val) in keys
      println(iob, "│   ", key, " = ", val)
    end
    println(iob, "└ @ ", something(_module, "nothing"), " ",
            something(filepath, "nothing"), ":", something(line, "nothing"))
  end
  write(logger.stream, take!(buf))
  nothing
end
