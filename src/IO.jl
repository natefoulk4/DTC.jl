
using Plots, FFTW

include("src.jl")

"
    plotter(effAvgResult; which={'worst', 'average', 'all'} )
Plot the autocorrelator of the DTC for the ``which`` qubit. ``'all'`` is not yet supported."
function plotter(effAvgResult; which="worst")
    allres, spinmap = effAvgResult
    worst_qubit = findmin(mean(allres[:,1:2:end,:], dims=(2,3))[:])[2]
    this_qubit = which == "worst" ? worst_qubit : which

    if typeof(this_qubit) == Int
        realres = allres[this_qubit,:]
    elseif which == "average"
        realres = mean(allres, dims=1)[:]
    elseif which == "all"
        println("FIXME. Plot all qubits.")
        return nothing
    end
    L = size(spinmap)[1]
    nperiods = size(spinmap)[2]
    res = abs.(fft(realres))
    if which == "all"
        println("FIXME. Plot all qubits")
        return nothing
    else
        plotly()
        l = @layout [a b ; c{0.2h}; d{0.2h}] 
        plot(realres, opacity=0.5, xscale=:log10, legend=false)
        plot!(collect(1:2:nperiods),realres[1:2:nperiods])
        p1 = plot!([1; collect(2:2:nperiods)],[-1.0 ; realres[2:2:nperiods]], title = string("Qubit #",this_qubit))
        p21 = plot(res, yscale =:log10, ylims=(0.1,350), legend=false)
        p3 = heatmap(spinmap[:,1:30], c=:viridis, clims=(-1,1))
        howitsgoing = nperiods < 1030 ? ((nperiods-29):nperiods) : 1001:1030
        p4 = heatmap(collect(howitsgoing), collect(1:L), spinmap[:,howitsgoing], c=:viridis, clims=(-1,1))#print(spinmap[:,1001:1030])
        plot(p1, p21, p3, p4, layout=l)
    end
end

plotter(effAvgResult1, effAvgResult2) = plotter((effAvgResult1, effAvgResult2))

"
    writeAutoCor(filename, output) 
writes the output of effAvgAutoCor (a 2-tuple) in binary"
function writeAutoCor(filename, output::Tuple)
    realres, spinmap = output
    open(filename, "w") do f
        # do realres
        write(f, Int64(length(realres)))
        for i in eachindex(realres)
            write(f, Float64(realres[i]))
        end

        # do spinmap
        write(f, Int64(size(spinmap)[1]))
        write(f, Int64(size(spinmap)[2]))
        for i in eachindex(spinmap)
            write(f, Float64(spinmap[i]))
        end

    end
    return nothing
end

"
    readAutoCor(filename) 
read the output of effAvgAutoCor (a 2-tuple) in binary"
function readAutoCor(filename)
    open(filename, "r") do f
        # do realres
        n = read(f, Int64)
        realres = zeros(Float64, n)
        for i in eachindex(realres)
            realres[i] = read(f, Float64)
        end

        # do spinmap
        m = read(f, Int64)
        n = read(f, Int64)
        spinmap = zeros(Float64, (m,n))
        for i in eachindex(spinmap)
            spinmap[i] = read(f, Float64)
        end

        return realres, spinmap
    end
end

"
    writeVector(filename, vector) 
Write a vector of floats in binary"
function writeVector(filename, vector)
    open(filename, "w") do f
        write(f, Int64(length(vector)))
        for i in eachindex(vector)
            write(f, Float64(vector[i]))
        end
    end
    return nothing
end

"
    readVector(filename) 
Read a vector of floats in binary"
function readVector(filename)
    open(filename, "r") do f

        n = read(f, Int64)
        vector = zeros(Float64, n)
        for i in eachindex(vector)
            vector[i] = read(f, Float64)
        end
        return vector
    end
end


"
    binToBase10(x::Vector{Int})
Convert a binary string into a base-10 number."
function binToBase10(x::Vector{Int64})
    L = length(x)
    basis = [2^(L-i) for i in 1:L]
    return(basis⋅x)
end