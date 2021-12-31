
using BenchmarkTools, Plots, FFTW

include("src.jl")

gr()
function plotter(realres, spinmap)
    nperiods = size(spinmap)[2]
    res = abs.(fft(realres))
    l = @layout [a b ; c{0.2h}; d{0.2h}] 
    plot(realres, opacity=0.5, xscale=:log10, legend=false)
    plot!(collect(1:2:nperiods),realres[1:2:nperiods])
    p1 = plot!([1; collect(2:2:nperiods)],[-1.0 ; realres[2:2:nperiods]])
    p21 = plot(res, yscale =:log10, ylims=(0.1,350), legend=false)
    p3 = heatmap(spinmap[:,1:30], c=:viridis, clims=(-1,1))
    howitsgoing = nperiods < 1030 ? ((nperiods-29):nperiods) : 1001:1030
    p4 = heatmap(collect(howitsgoing), collect(1:length(init)), spinmap[:,howitsgoing], c=:viridis, clims=(-1,1))#print(spinmap[:,1001:1030])
    plot(p1,p21, p3, p4, layout=l)
end

plotter(x::Tuple) = plotter(x[1], x[2])

function savedata(filename, output::Tuple)
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

function readdata(filename)
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