using Unitful
using Unitful: ms, s, kHz, Hz

TimeDim = Unitful.Dimensions{(Unitful.Dimension{:Time}(1//1),)}
SecondsUnit = Unitful.FreeUnits{(Unitful.Unit{:Second,TimeDim}(0,1//1),),
                                      TimeDim}
Seconds{N} = Quantity{N,TimeDim,SecondsUnit}
