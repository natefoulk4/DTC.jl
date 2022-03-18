module DTC

using LinearAlgebra, Statistics, Distributions, SparseArrays, TimerOutputs, Plots, FFTW

export getIsingNNJtensor, autocorrelator, effAvgAutoCor, plotter, LsrsOverParamRange

const to = TimerOutput()

include("IO.jl")
include("src.jl")

end # module
