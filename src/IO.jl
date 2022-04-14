
using Plots, FFTW


"
    plotter(effAvgResult; which={Int, 'worst', 'average', 'all'} )
Plot the autocorrelator of the DTC for the ``which`` qubit. ``'all'`` is not yet supported."
function plotter(effAvgResult; which="worst")
    if length(size(effAvgResult[1])) == 2
        allres, spinmap = effAvgResult
        worst_qubit = findmin(mean(allres[:,1:2:end,:], dims=(2,3))[:])[2]
        this_qubit = which == "worst" ? worst_qubit : which

        if typeof(this_qubit) == Int
            realres = allres[this_qubit,:]
        elseif lowercase(which) == "average"
            realres = mean(allres, dims=1)[:]
        elseif lowercase(which) == "all"
            println("FIXME. Plot all qubits.")
            return nothing
        end
    elseif length(size(effAvgResult[1])) == 1 # Legacy support
        realres, spinmap = effAvgResult
        L = size(spinmap)[1]
        this_qubit = Int(round(L/2))
    end
    L = size(spinmap)[1]
    nperiods = size(spinmap)[2]
    res = abs.(fft(realres))
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

function writeArray(filename, arr)
    open(filename, "w") do f
        d = length(size(arr))
        write(f, Int64(d))
        for j in 1:d
            write(f, Int64(size(arr)[d]))
        end
        for i in eachindex(arr)
            write(f, Float64(arr[i]))
        end
    end
    return nothing
end

function readArray(filename)
    open(filename, "r") do f
        d = read(f, Int64)
        sizeVec = zeros(Int64, d)
        for j in 1:d
            sizeVec[j] = read(f, Int64)
        end
        sizeTup = Tuple(sizeVec)
        arr = zeros(sizeTup)
        for i in eachindex(arr)
            arr[i] = read(f, Float64)
        end
        return arr
    end
end

