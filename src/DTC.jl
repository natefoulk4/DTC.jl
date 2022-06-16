module DTC

using LinearAlgebra, Statistics, Distributions, SparseArrays, TimerOutputs, Plots, FFTW, LaTeXStrings

export autocorrelator, effAvgAutoCor, plotter, LsrsOverParamRange, readArray, writeArray, readAutoCor, writeAutoCor, binToBase10, base10ToBin, levelspacing, avgLevelSpacings, U1, matrix_density, to, MAX_PERIODS, logRange, logIntRange

export getFig1Data, getFig2Data, getFig3Data, getFig4Data, getFig5Data, getFig6Data, getFig7Data, Figure1and2, Figure3and4, Figure5, Figure6

const to = TimerOutput()

include("IO.jl")
include("src.jl")
include("makeFigures.jl")

end # module
