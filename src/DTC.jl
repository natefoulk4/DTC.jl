module DTC

using LinearAlgebra, Statistics, Distributions, SparseArrays, TimerOutputs, Plots, FFTW

export autocorrelator, effAvgAutoCor, plotter, LsrsOverParamRange, readArray, writeArray, readAutoCor, writeAutoCor, binToBase10, base10ToBin, levelspacing, avgLevelSpacings, U1, matrix_density, to

const to = TimerOutput()

include("IO.jl")
include("src.jl")

end # module
