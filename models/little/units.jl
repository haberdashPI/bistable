using Unitful
using Unitful: ms, s, kHz, Hz

TimeDim = Unitful.Dimensions{(Unitful.Dimension{:Time}(1//1),)}
SecondsUnit = Unitful.FreeUnits{(Unitful.Unit{:Second,TimeDim}(0,1//1),),TimeDim}
Seconds{N} = Quantity{N,TimeDim,SecondsUnit}

FreqDim = Unitful.Dimensions{(Unitful.Dimension{:Time}(-1//1),)}
HzUnit = Unitful.FreeUnits{(Unitful.Unit{:Hertz,FreqDim}(0,1//1),),FreqDim}
Hertz{N} = Quantity{N,FreqDim,HzUnit}
