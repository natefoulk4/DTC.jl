module DTC

using LinearAlgebra, Statistics, Distributions, SparseArrays, TimerOutputs, Plots, FFTW

export getIsingNNJtensor, σx

const to = TimerOutput()

include("IO.jl")
include("src.jl")

end # module
