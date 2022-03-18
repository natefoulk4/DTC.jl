
using LinearAlgebra, Statistics, Distributions, TimerOutputs

"
    Rx(n, φ, spinKet, newKet)
Operate on ``spinKet`` (length ``2^L``) with the rotation unitary R^x_n(φ) (on the ``n``th qubit). Return ``newKet``"
@timeit to function Rx(n::Int, φ::Real, spinKet::Vector{ComplexF64}, newKet::Vector{ComplexF64})
    L = Int(log2(length(spinKet)))
    stride =  2^(L-n)   # 1 if n=4, 2 if n=3, 4 if n=2, 8 if n = 1
    stride2 = 2^(n-1)   # 8 if n=4, 4 if n=3, 2 if n=2, 1 if n = 1

    cos1 = cos(φ)
    sin1 = sin(φ)

    for i in 1:stride2
        for j in 1:stride
            #println(2*(i-1)*stride+j, 2*(i-1)*stride+j + stride)
            newKet[2*(i-1)*stride+j] = cos1*spinKet[2*(i-1)*stride+j] + im*sin1*spinKet[2*(i-1)*stride+j + stride]
            newKet[2*(i-1)*stride+j + stride]  = im*sin1*spinKet[2*(i-1)*stride+j] + cos1*spinKet[2*(i-1)*stride+j + stride]
        end
    end

    return newKet
end

"
    getIsingNNJtensor(L)
Return ``2^L x L`` matrix, where the ``n``th column contains the approriate signs for the diagonal of the J\\_n matrix. _Ising model only._"
@timeit to function getIsingNNJtensor(L)
    diagonals = fill(1.0 + 0.0im, (2^L,L))
    for n in 1:L # n is the nth exchange coupling (J_{n, n+1})
        stride =  2^(L-n)   # 1 if n=4, 2 if n=3, 4 if n = 2, 8 if n=1
        halfstride = Int(round(stride/2))
        stride2 = 2^(n-1)   # 8 if n=4, 4 if n=3, 2 if n=2, 1 if n = 1
        # L = 2
        # dd du ud uu
        # 1 -1 -1 1    N = 1
        # 1 -1 -1 1    N = 2
        # L = 3 
        # ddd ddu dud duu udd udu uud uuu
        # 1  1 -1 -1 -1 -1  1  1   4,2,1   N=1 [5-12]
        # 1 -1 -1  1  1 -1 -1  1   2,1,2   N=2 [3-6], [11-14]
        # 1 -1  1 -1 -1  1 -1  1   1,1,4   N=3 [2-3], [6-7], [10-11], [14-15]
        # L = 4
        # dddd dddu ddud dduu dudd dudu duud duuu uddd ....
        # 1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1   8,4,1   N=1 [5-12]
        # 1  1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1  1  1   4,2,2   N=2 [5-12]
        # 1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1   2,1,4   N=3 [3-6], [11-14]
        # 1 -1  1 -1  1 -1  1 -1 -1  1 -1  1 -1  1 -1  1   1,1,8   N=4 [2-3], [6-7], [10-11], [14-15]
        #
        for i in 1:stride2
            for j in 1:stride
                diagonals[2*(i-1)*stride + halfstride + j, n] *= -1.0
            end
        end
    end

    #Fix the n=L case
    diagonals[1:2^(L-1), L] *= -1.0

    return diagonals
end

"
    getHtensor(L)
Return a ``2^L x L`` matrix where the ``n``th column contains the appropriate signs depending on whether that ``n``th spin is up or down for that base ket."
function getHtensor(L)
    hTensor = fill(1.0+0.0im, (2^L, L))
    for n in 1:L    
        stride =  2^(L-n)   # 1 if n=4, 2 if n=3, 4 if n=2, 8 if n = 1
        stride2 = 2^(n-1)   # 8 if n=4, 4 if n=3, 2 if n=2, 1 if n = 1
    
        for i in 1:stride2
            for j in 1:stride
               hTensor[2*(i-1)*stride+j, n] = -1.0+0.0im
            end
        end
    end
    return hTensor
end

"
    levelspacing(vals)
Calculate level spacing ratios (LSRs) of a list of eigenvalues (not necessarily sorted). Return mean of the LSRs."
@timeit to function levelspacing(vals)
    sort!(vals)
    for i in 1:length(vals)-1
        vals[i] = vals[i+1] - vals[i]
    end
    for i in 1:length(vals)-2
        vals[i] = min(vals[i],vals[i+1])/max(vals[i],vals[i+1])
    end
    return mean(vals[1:end-2])
end

"
    efficientHam(Hspace, hs, js, jMat, hTensor; BCs)
Multiply each J and H columns of ``jMat`` and ``hTensor`` by the appropriate prefactors contained in ``js`` and ``hs``. Add them together and assign them to the diagonal of ``Hspace``, overwriting ``Hspace`` completely. Return ``Hspace``. keyword: ``BCs={'open', 'periodic'}``"
@timeit to function efficientHam(Hspace, hs, js, jMat, hTensor; BCs)
    L = Int(log2(size(Hspace)[1]))
    for i in 1:2^L
        Hspace[i,i] = 0.0 + 0.0im
    end
    if BCs == "open"
        for j in 1:L-1, n in 1:2^L
            Hspace[n,n] += (js[j] * jMat[n,j])
        end
    elseif BCs == "periodic"
        for j in 1:L, n in 1:2^L
            Hspace[n,n] += (js[j] * jMat[n,j])
        end
    end
    for j in 1:L, n in 1:2^L
        Hspace[n,n] += hs[j] * hTensor[n,j]
    end
    return Hspace
end

"
    getKet(spinArray)
Given a (length ``L``) array of spins (must be a pure state), find the appropriate base ket (length ``2^L``) of Hilbert space."
function getKet(spinArray)
    L = length(spinArray)
    spinBasis=[reverse(digits(i, base=2, pad=length(spinArray))) for i in 0:2^length(spinArray)-1]
    coeffs = zeros(ComplexF64,2^L)
    for i in eachindex(spinBasis)
        if spinBasis[i] == spinArray  
            coeffs[i] = 1.0 + 0.0im
        end
    end
    return coeffs
end

"
    getSpins!(ket, spinBasis, spinMatrix, ind)
Convert the base ``ket``` (length ``2^L``) to an array of spins (projected to z-axis). Write these spins in the ``ind``th column of ``spinMatrix``."
@timeit to function getSpins!(ket, spinBasis, spinMatrix, ind)
   for i in 1:size(spinBasis)[1] #2^L
        for j in 1:size(spinBasis)[2] #L
            spinMatrix[j,ind] += abs2(ket[i])*spinBasis[i,j]
        end
    end
    return nothing
end

"
    IsingefficU2(Hspace, hs, jz, jArrays, hTensor; BCs='open')
Given memory ``Hspace``, call ``efficientHam`` and then exponentiate the diagonal elements in-place. Return Hspace as U2, the MBL unitary."
@timeit to function IsingefficU2(Hspace, hs, js, jMat, hTensor; BCs="open") 
    efficientHam(Hspace, hs, js, jMat, hTensor; BCs=BCs)
    for i in 1:size(Hspace)[1]
        Hspace[i,i] = exp(-im*Hspace[i,i])
    end
    return Hspace
    #return exp(-im.*efficientHam(Hspace, hs, js, jTensor, hTensor))
end

"
    matrix_density(mat)
Given a matrix ``mat``, return the density of the matrix, which is the ratio of nonzero entries to the total number of matrix entries."
matrix_density(mat::Matrix) = length(findall(!iszero,mat))/length(mat)

"
    autocorrelator(spins, basis, eps, U2, N)
Simulate the dynamics of a DTC for a given instance of disorder. Take initial state of ``spins`` (length ``L``) and ``basis`` (size ``2^L x L``). Apply Floquet unitaries (using ``eps``, the perturbation of a π pulse, and ``U2``) ``N`` times. Return two matrices both size ``L x N+1``. 1st return (autoCor) is the correlation relative to the initial state. The second (moreSpins) is the _absolute_ spin of each qubit"
@timeit to function autocorrelator(spins, basis, eps, U2, N)
    initKet = getKet(spins)
    L = length(spins)
    negOneSpins = replace(spins, 0 => -1)


    autoCor = zeros(L,N+1)
    moreSpins = zeros(L,N+1)
    currentKet = deepcopy(initKet)
    interKet = zeros(ComplexF64, 2^L)
    newKet = zeros(ComplexF64, 2^L)

    autoCor[:,1] .= 1.0
    moreSpins[:,1] = negOneSpins
     for i in 2:N+1
        for k in 1:L
            #mul!(interKet, U1s[k], currentKet)
            Rx(k, pi/2*(1-eps), currentKet, interKet)
            interKet, currentKet = currentKet, interKet
        end

        @timeit to "U2 mul" mul!(newKet,U2,currentKet)
    
        getSpins!(newKet, basis, moreSpins, i)

        currentKet, newKet = newKet, currentKet
    end

    autoCor = moreSpins .* negOneSpins

    return autoCor, moreSpins
end

"
    getHsAndJs(niters, L, J0, σJ, σH)
Sample distribution ``niter*L`` times. Return Hmatrix (``niters x L``) and J matrix (``niters x L``)"
function getHsAndJs(niters, L, J0, σj, σh)
    hs = zeros(niters, L)
    js = zeros(niters, L)
    M = 8
    discrete_Jrange = J0 .+ (σj .* cos.(pi*collect(0:M-1)/(M-1)))

    for i in eachindex(hs)
        if σh > 0.0
            hs[i] = rand(Normal(0.0, σh))
        else
            hs[i] = 0.0
        end
    end

    for i in eachindex(js)
        if σj > 0.0
            js[i] = rand(discrete_Jrange)
        else
            js[i] = J0
        end
    end

    return hs, js
end

"
    effAvgAutoCor(niters, nperiods, spins, ε, J0, σJ, σH; t, BCs)
Exact same as autocorrelator, except average over ``niters`` simulations. Return ``finalCors`` and ``allSpins``, which are both ``2^L x L`` matrices "
@timeit to function effAvgAutoCor(niters, nperiods, spins, ε, J0, σj, σh; t=0.0, BCs="open")
    L = length(spins)
    Hspace = spzeros(ComplexF64, 2^L, 2^L)

    hs, js = getHsAndJs(niters, L, J0, σj, σh)

    
    cors = zeros(L, nperiods+1, niters)
    finalCors=zeros(L, nperiods+1)
    allSpins = zeros(L, nperiods+1, niters)
    jIsingTensor = getIsingNNJtensor(L)
    hTensor = getHtensor(L)            # more of the time

    basis = zeros(2^L, L)
    for i in 1:2^L
        basis[i,:] = Float64.(reverse(digits(i-1, base=2, pad=L)))
    end
    replace!(x->iszero(x) ? -1.0 : x, basis) #This basis is a matrix of the spins
    
    #u1s = betterU1(L, ε)

    for i in 1:niters
        cors[:,:,i],allSpins[:,:, i]  = autocorrelator(spins, basis, ε, IsingefficU2(Hspace,  hs[i,:] ,  js[i,:], jIsingTensor, hTensor; BCs=BCs), nperiods)
        if i % (niters/10) == 0
            println("Finished ",i,"th iteration")
        end
    end

    finalCors = mean(cors, dims=3)
    allSpins[:,:,1] = mean(allSpins,dims=3)

    return (finalCors, allSpins[:,:,1])
end

"
    avgLevelSpacings(niters, nperiods, spins, ε, J0, σJ, σH; t=0.0, BCs='open')
Calculate ``niters`` of the level spacing ratios. Return the average LSR overall."
@timeit to function avgLevelSpacings(niters, nperiods, spins, ε, J0, σj, σh; t=0.0, BCs="open")
    L = length(spins)
    Hspace = zeros(ComplexF64, 2^L, 2^L)
    Hspace2 = Hermitian(zeros(ComplexF64, 2^L, 2^L))

    hs, js = getHsAndJs(niters, L, J0, σj, σh)
    rats = zeros(niters)
    jIsingTensor = getIsingNNJtensor(L)
    hTensor = getHtensor(L)            # more of the time

    basis = zeros(2^L, L)
    for i in 1:2^L
        basis[i,:] = Float64.(reverse(digits(i-1, base=2, pad=L)))
    end
    replace!(x->iszero(x) ? -1.0 : x, basis) #This basis is a matrix of the spins
    
    u1 = newU1(L, ε) #FIXME
    noconverge=0

    for i in 1:niters
        @timeit to "matmul" Hspace2 = Array(u1*IsingefficU2(Hspace,hs[i,:],js[i,:],jIsingTensor,hTensor, BCs=BCs))
        @timeit to "eigs" vals = eigvals!(Hspace2)
        #vals = Real.(diag(Array(efficientHam(Hspace,hs[i,:],js[i,:],jIsingTensor,hTensor; BCs=BCs))))
        #println(vals)
        rats[i] = levelspacing( mod.(Real.(round.(log.(vals) .* im,digits=8)), 2*pi) )
        #rats[i] = levelspacing( vals )
        if i % (niters/10) == 0
            #println("Finished $i th iteration out of $niters")
        end
    end

    return mean(rats)
end

"
    LsrsOverParamRange(param, Lrange; BCs)
Calculate LSRs over a parameter range (``param={'eps','j','sigH'}``). Done for several values of L (``Lrange``). Boundary conditions (``BCs={'open','periodic'}``) optional."
function LsrsOverParamRange(param, Lrange; BCs)

    if param == "eps"
        paramRange = range(0.0, 1.0, step=0.05)
    elseif param == "j"
        paramRange = 10 .^ collect(range(-2,1.5;step=0.2))
    elseif param == "sigH"
        paramRange = 10 .^ collect(range(-2,2;step=0.4))
    else
        error("Parameter not correctly specificied!")
    end
    n = length(paramRange)
    lsrs = zeros(length(Lrange), n)

    niter = 1000
    nperiods = 1000
    J0 = pi/4
    σJ = pi/2
    σH  = pi/50
    ε = 0.1
    i = 1

    for l in Lrange
        println("Starting L=",l)
        for (j,param) in enumerate(paramRange) 
            if param == "eps"
                ε = param
            elseif param == "j"
                J0 = param
            elseif param == "sigH"
                σH = param
            end
            init = rand([0,1], l)
            lsrs[i,j] = avgLevelSpacings(niter, nperiods, init, ε, J0, σJ, σH ; BCs=BCs)
            println("Finished $j th data point out of $n ($niter iterations each point)")
        end
        i += 1
    end

    return paramRange, lsrs'
end