using DTC, Statistics, LaTeXStrings, Plots, Plots.PlotMeasures
#pyplot()
 
function getFig1Data()
    xs = zeros(20)
    ys = zeros(20)
    data = zeros(20,20,8)
    εs, js = range(0.0, 0.33, length=20), range(0.0, pi, length=20)

    inits = [[1,0,0,0]];
    ns = [1,3];
    nperiods = [2, 10, 50, 200];
    ind = 0
    for (i,init) in enumerate(inits), n in ns, nperiod in nperiods
        ind = ind + 1
        for (e, thisEps) in enumerate(εs), (j,thisJ0) in enumerate(js)
            xs[j] = thisJ0
            ys[e] = thisEps
            data[e,j,ind] = minimum(abs.(effAvgAutoCor(100, nperiod, init; ε=thisEps, J0=thisJ0, σj=0.0, h0=20000, σh=50, H_dist="uniform", J_dist="uniform", t=0.0, num_H2I=0, BCs="open", diagonalization=false)[1][n,1:2:end]))
        end
    end
    return xs, ys, data
end
function getFig2Data()
    xs = zeros(20)
    ys = zeros(20)
    data = zeros(20,20,8)
    εs, js = range(0.0, 0.33, length=20), range(0.0, pi, length=20)

    inits = [[1,0,0,0]];
    ns = [1,3];
    nperiods = [2, 10, 50, 200];
    ind = 0
    for (i,init) in enumerate(inits), n in ns, nperiod in nperiods
        ind = ind + 1
        for (e, thisEps) in enumerate(εs), (j,thisJ0) in enumerate(js)
            xs[j] = thisJ0
            ys[e] = thisEps
            data[e,j,ind] = minimum(abs.(effAvgAutoCor(100, nperiod, init; ε=thisEps, J0=thisJ0, σj=3.0, h0=20000, σh=50, H_dist="uniform", J_dist="uniform", t=0.0, num_H2I=0, BCs="open", diagonalization=false)[1][n,1:2:end]))
        end
    end
    return xs, ys, data
end
function getFig3Data()
    xs = zeros(20)
    ys = zeros(20)
    data = zeros(20,20,8)
    εs, js = range(0.0, 0.33, length=20), range(0.0, pi, length=20)

    inits = [[1,0,1,0], [1,0,0,0], [0,0,1,0], [0,1,1,0], [0,0,0,0], [1,1,0,0], [1,1,1,1], [1,1,1,0]];
    ns = [3];
    nperiods = [200];
    ind = 0
    for (i,init) in enumerate(inits), n in ns, nperiod in nperiods
        ind = ind + 1
        for (e, thisEps) in enumerate(εs), (j,thisJ0) in enumerate(js)
            xs[j] = thisJ0
            ys[e] = thisEps
            data[e,j,ind] = minimum(abs.(effAvgAutoCor(100, nperiod, init; ε=thisEps, J0=thisJ0, σj=0.0, h0=20000, σh=50, H_dist="uniform", J_dist="uniform", t=0.0, num_H2I=0, BCs="open", diagonalization=false)[1][n,1:2:end]))
        end
    end
    return xs, ys, data
end
function getFig4Data()
    xs = zeros(20)
    ys = zeros(20)
    data = zeros(20,20,8)
    εs, js = range(0.0, 0.33, length=20), range(0.0, pi, length=20)

    inits = [[1,0,1,0], [1,0,0,0], [0,0,1,0], [0,1,1,0], [0,0,0,0], [1,1,0,0], [1,1,1,1], [1,1,1,0]];
    ns = [3];
    nperiods = [200]; 
    ind = 0
    for (i,init) in enumerate(inits), n in ns, nperiod in nperiods
        ind = ind + 1
        for (e, thisEps) in enumerate(εs), (j,thisJ0) in enumerate(js)
            xs[j] = thisJ0
            ys[e] = thisEps
            data[e,j,ind] = minimum(abs.(effAvgAutoCor(100, nperiod, init; ε=thisEps, J0=thisJ0, σj=3.0, h0=20000, σh=50, H_dist="uniform", J_dist="uniform", t=0.0, num_H2I=0, BCs="open", diagonalization=false)[1][n,1:2:end]))
        end
    end
    return xs, ys, data
end
function getFig5Data()
    xNum = 10
    yNum = 10
    xs = zeros(xNum)
    ys = zeros(yNum)
    data = zeros(xNum,yNum,2)
    εs, σjs = range(0.0, 0.10, length=xNum), 10 .^ range(-2.5, 2.0, length=yNum)

    inits = [[1,0,0,0]];
    ns = [3];
    nperiods = [2000]; 
    j0s = [1.5, 10000];
    ind = 0
    for (i,init) in enumerate(inits), n in ns, nperiod in nperiods, j0 in j0s
        ind = ind + 1
        for (e, thisEps) in enumerate(εs), (j,thisJ0) in enumerate(σjs)
            xs[j] = thisJ0
            ys[e] = thisEps
            data[e,j,ind] = minimum(abs.(effAvgAutoCor(1000, nperiod, init; ε=thisEps, J0=j0, σj=thisJ0, h0=20000, σh=50, H_dist="uniform", J_dist="uniform", t=0.0, num_H2I=0, BCs="open", diagonalization=false)[1][n,1:2:end]))
        end
    end
    return xs, ys, data
end
function getFig6Data()
    xNum = 10
    yNum = 10
    xs = zeros(xNum)
    ys = zeros(yNum)
    data = zeros(xNum,yNum,2)
    εs, σhs = range(0.0, 0.10, length=xNum), 10 .^ range(-2.5, 2.0, length=yNum)

    inits = [[1,0,0,0]];
    ns = [3];
    nperiods = [2000]; 
    h0s = [1.5, 10000];
    ind = 0
    for (i,init) in enumerate(inits), n in ns, nperiod in nperiods, h0 in h0s
        ind = ind + 1
        for (e, thisEps) in enumerate(εs), (j,thisJ0) in enumerate(σhs)
            xs[j] = thisJ0
            ys[e] = thisEps
            data[e,j,ind] = minimum(abs.(effAvgAutoCor(1000, nperiod, init; ε=thisEps, J0=1.5, σj=3, h0=h0, σh=thisJ0, H_dist="uniform", J_dist="uniform", t=0.0, num_H2I=0, BCs="open", diagonalization=false)[1][n,1:2:end]))
        end
    end
    return xs, ys, data
end

function Figure1and2(x,y,z)
    myPlots = Array{Any}(undef, 8)
    fo = ("times", 20)
    Ts = [2, 10, 50, 200, 2, 10, 50 , 200]
    for i in 1:8
        yt = xt = []
        yl = xl = ""
        bm = 0mm
        lm = 5mm
        if i == 1 || i == 5
            yl = L"\varepsilon"
            yt = [0,0.1,0.2,0.3]
            lm = 8mm
        end
        if i > 4
            xl = L"J_0"
            xt = ([0, pi/4, pi/2, 3*pi/4, pi], [L"0", L"\frac{\pi}{4}", L"\frac{\pi}{2}", L"\frac{3\pi}{4}", L"\pi"])
            bm = 12mm
        end
        myPlots[i] = contourf(x,y,z[:,:,i], clim=(0,1), colorbar=false, title=latexstring(" t = ",string(Ts[i]...),"T"), c=:thermal, levels=9, yticks=yt, xticks=xt, ylabel=yl, xlabel=xl, tickfont=("times",15), guidefont=("times", 25), titlefont=fo, tick_direction=:out, topmargin=4mm, bottommargin=bm, leftmargin=lm)
    end


    h1 = scatter([0,0], [0,1], zcolor=[0,3], clims=(0,1), label="", c=palette(:thermal,10), colorbar_title=L"Z_n(t)", colorbar_titlefont=("serif-roman", 20), tickfont=("times", 15), xlims=(1.0,1.1), framestyle=:none, axis=false, right_margin=5mm, left_margin=5mm)

    h2 = scatter([0,0], [0,1], annotations=[(1.04, 0.78, (L"n=1", 20)),(1.04, 0.22, (L"n=3", 20))],  zcolor=[0,3], clims=(0,1), label="", c=palette(:thermal,10), colorbar_title="", colorbar=false,  xlims=(1.0,1.1), axis=false, framestyle=:none, ticks=[])

    l = @layout [ (2,4) a{0.06w} b{0.04w}]

    return plot(myPlots..., h2, h1, layout=l, size=(1500,750))
end
function Figure3and4(x,y,z)
    myPlots = Array{Any}(undef, 8)
    fo = ("times", 20)
    inits = [[1,0,1,0], [1,0,0,0], [0,0,1,0], [0,1,1,0], [0,0,0,0], [1,1,0,0], [1,1,1,1], [1,1,1,0]];
    for i in 1:8
        yt = xt = []
        yl = xl = ""
        bm = 0mm
        lm = 5mm
        if i == 1 || i == 5
            yl = L"\varepsilon"
            yt = [0,0.1,0.2,0.3]
            lm = 8mm
        end
        if i > 4
            xl = L"J_0"
            xt = ([0, pi/4, pi/2, 3*pi/4, pi], [L"0", L"\frac{\pi}{4}", L"\frac{\pi}{2}", L"\frac{3\pi}{4}", L"\pi"])
            bm = 12mm
        end
        myPlots[i] = contourf(x,y,z[:,:,i], clim=(0,1), colorbar=false, title=latexstring(" | \\psi_0 \\rangle = | ",string(inits[i]...), " \\rangle"), c=:thermal, levels=9, yticks=yt, xticks=xt, ylabel=yl, xlabel=xl, tickfont=("times",15), guidefont=("times", 25), titlefont=fo, tick_direction=:out, topmargin=4mm, bottommargin=bm, leftmargin=lm)
    end


    h1 = scatter([0,0], [0,1], zcolor=[0,3], clims=(0,1), label="", c=palette(:thermal,10), colorbar_title=L"Z_3(t)", colorbar_titlefont=("serif-roman", 20), tickfont=("times", 15), xlims=(1.0,1.1), framestyle=:none, axis=false, right_margin=5mm, left_margin=5mm)

    l = @layout [ (2,4) b{0.04w}]

    return plot(myPlots..., h1, layout=l, size=(1500,750))
end
function Figure5(x,y,z)
    myPlots = Array{Any}(undef, 2)
    fo = ("serif", 25)
    yl = L"\varepsilon"
    yt = [0,0.03,0.06,0.1]
    xl = L"\sigma_J"
    xt = [10^(-2), 10^(-1), 10^(0), 10^(1), 10^(2)] #, [L"0", L"\frac{\pi}{4}", L"\frac{\pi}{2}", L"\frac{3\pi}{4}", L"\pi"])

    myPlots[1] = contourf(x,y,z[:,:,1], colorbar=false, c=:thermal, levels=9, yticks=yt, xticks=xt, ylabel=yl, xlabel=xl, tickfont=("serif-roman",10), guidefont=fo, titlefont=fo, xscale=:log10, topmargin=5mm, left_margin=6mm, bottommargin=10mm, title=L" J_0 = 1.5 ")
    myPlots[2] = contourf(x,y,z[:,:,2], colorbar=false, c=:thermal, levels=9, xticks=xt, yticks=yt, ylabel=yl, xlabel=xl, tickfont=("serif-roman",10), guidefont=fo, titlefont=fo, xscale=:log10, topmargin=5mm, left_margin=6mm, bottommargin=10mm, title=L" J_0 = 10,\!000 ")

    h1 = scatter([0,0], [0,1], zcolor=[0,3], clims=(0,1), label="", c=palette(:thermal,10), colorbar_title=L"Z_3(t)", colorbar_titlefont=("serif-roman", 20), tickfont=("times", 12), xlims=(1.0,1.1), framestyle=:none, axis=false, right_margin=5mm, left_margin=5mm)

    l = @layout [Plots.grid(1,2) a{0.04w}]

    return plot(myPlots..., h1, layout=l, size=(1000,500))
end
function Figure6(x,y,z)
    myPlots = Array{Any}(undef, 2)
    fo = ("serif", 25)
    yl = L"\varepsilon"
    yt = [0,0.03,0.06,0.1]
    xl = L"\sigma_h"
    xt = [10^(-2), 10^(-1), 10^(0), 10^(1), 10^(2)]#, [L"10^{-2}", L"10^{-1}", L"10^0", L"10^{1}", L"10^{2}"])

    myPlots[1] = contourf(x,y,z[:,:,1], colorbar=false, c=:thermal, levels=9, yticks=yt, xticks=xt, ylabel=yl, xlabel=xl, tickfont=("serif-roman",10), guidefont=fo, titlefont=fo, xscale=:log10, topmargin=5mm, left_margin=6mm, bottommargin=10mm, title=L" h_0 = 1.5")
    myPlots[2] = contourf(x,y,z[:,:,2], colorbar=false, c=:thermal, levels=9, xticks=xt, yticks=yt, ylabel=yl, xlabel=xl, tickfont=("serif-roman",10), guidefont=fo, titlefont=fo, xscale=:log10, topmargin=5mm, left_margin=6mm, bottommargin=10mm, title=L" h_0 = 10,\!000")

    h1 = scatter([0,0], [0,1], zcolor=[0,3], clims=(0,1), label="", c=palette(:thermal,10), colorbar_title=L"Z_3(t)", colorbar_titlefont=("serif-roman", 20), tickfont=("times", 12), xlims=(1.0,1.1), framestyle=:none, axis=false, right_margin=5mm, left_margin=5mm)

    l = @layout [a b c{0.04w}]

    return plot(myPlots..., h1, layout=l, size=(1000,500))
end