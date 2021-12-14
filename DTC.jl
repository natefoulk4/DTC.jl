
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