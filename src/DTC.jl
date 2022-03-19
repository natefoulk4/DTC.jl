module DTC

using LinearAlgebra, Statistics, Distributions, SparseArrays, TimerOutputs, Plots, FFTW

export getIsingNNJtensor, autocorrelator, effAvgAutoCor, plotter, LsrsOverParamRange, readVector, writeVector, readAutoCor, writeAutoCor, binToBase10, base10ToBin, levelspacing, avgLevelSpacings

const to = TimerOutput()

include("IO.jl")
include("src.jl")

end # module
