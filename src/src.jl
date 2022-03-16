
using LinearAlgebra, Statistics, Distributions, SparseArrays, TimerOutputs

function initialize(spins::Vector{Vector{Int64}}, coeffs::SparseVector{ComplexF64, Int64})
    L = length(spins[1])
    zeroSpinor = ( [(reverse(digits(i, base=2, pad=L))) for i in 0:2^L-1], spzeros(ComplexF64, 2^L) )
    for i in 1:length(spins)
        for j in 1:length(zeroSpinor[1])
            if spins[i] == zeroSpinor[1][j]
                zeroSpinor[2][j] += coeffs[i]
            end
        end
    end
    return zeroSpinor[2]
end
function σx(n::Int, spinArray::SparseVector{ComplexF64, Int64}, newArray::SparseVector{ComplexF64, Int64})
    L = Int(log2(length(spinArray)))
    stride =  2^(L-n)   # 1 if n=4, 2 if n=3, 4 if n=2, 8 if n = 1
    stride2 = 2^(n-1)   # 8 if n=4, 4 if n=3, 2 if n=2, 1 if n = 1

    for i in 1:stride2
        for j in 1:stride
           newArray[2*(i-1)*stride+j], newArray[2*(i-1)*stride+j + stride] = spinArray[2*(i-1)*stride+j + stride], spinArray[2*(i-1)*stride+j]
        end
    end

    return newArray
end
@timeit to function Rx(n::Int, θ::Real, spinArray::Vector{ComplexF64}, newArray::Vector{ComplexF64})
    L = Int(log2(length(spinArray)))
    stride =  2^(L-n)   # 1 if n=4, 2 if n=3, 4 if n=2, 8 if n = 1
    stride2 = 2^(n-1)   # 8 if n=4, 4 if n=3, 2 if n=2, 1 if n = 1

    cos1 = cos(θ)
    sin1 = sin(θ)

    for i in 1:stride2
        for j in 1:stride
            #println(2*(i-1)*stride+j, 2*(i-1)*stride+j + stride)
            newArray[2*(i-1)*stride+j] = cos1*spinArray[2*(i-1)*stride+j] + im*sin1*spinArray[2*(i-1)*stride+j + stride]
            newArray[2*(i-1)*stride+j + stride]  = im*sin1*spinArray[2*(i-1)*stride+j] + cos1*spinArray[2*(i-1)*stride+j + stride]
        end
    end

    return newArray
end

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

function σy(n::Int, spinArray::SparseVector{ComplexF64, Int64}, newArray::SparseVector{ComplexF64, Int64})
    L = Int(log2(length(spinArray)))
    stride =  2^(L-n)   # 1 if n=4, 2 if n=3, 4 if n=2, 8 if n = 1
    stride2 = 2^(n-1)   # 8 if n=4, 4 if n=3, 2 if n=2, 1 if n = 1

    for i in 1:stride2
        for j in 1:stride
           newArray[2*(i-1)*stride+j], newArray[2*(i-1)*stride+j + stride] = im*spinArray[2*(i-1)*stride+j + stride], -im*spinArray[2*(i-1)*stride+j]
        end
    end

    return newArray
end

function σz(n::Int, spinArray::SparseVector{ComplexF64, Int64}, newArray::SparseVector{ComplexF64, Int64})
    L = Int(log2(length(spinArray)))
    stride =  2^(L-n)   # 1 if n=4, 2 if n=3, 4 if n=2, 8 if n = 1
    stride2 = 2^(n-1)   # 8 if n=4, 4 if n=3, 2 if n=2, 1 if n = 1

    for i in 1:stride2
        for j in 1:stride
           newArray[2*(i-1)*stride+j], newArray[2*(i-1)*stride+j + stride] = -spinArray[2*(i-1)*stride+j], spinArray[2*(i-1)*stride+j + stride]
        end
    end

    return newArray
end

function efficσiσj(i,j, spinor1, spinor2, spinor3, theta) # spinor2 and 3 just need to be the same size
    
    if theta != 0.0
        σx(i,σx(j,spinor1, spinor2),spinor2)
        for i in 1:length(spinor3) 
            spinor3[i] = theta*spinor2[i]
        end
        
        σy(i,σy(j,spinor1, spinor2),spinor2)
        for i in 1:length(spinor3) 
            spinor3[i] += theta*spinor2[i]
        end

        σz(i,σz(j,spinor1, spinor2),spinor2)
        for i in 1:length(spinor3) 
            spinor3[i] += spinor2[i]
        end
    else
        σz(i,σz(j,spinor1, spinor2),spinor2)
        for i in 1:length(spinor3) 
            spinor3[i] = spinor2[i]
        end
    end
    return spinor3
    #return theta*σx!(i,σx!(j, spinor)) + theta*σy!(i,σy!(j, spinor)) + σz!(i,σz!(j, spinor))
end

function getBasis(L) 
    bas = [spzeros(ComplexF64,2^L) for i in 1:2^L]
    for i in 1:2^L
        bas[i][i] = 1.0 + 0.0im
    end
    return bas
end

function operatorToMatrix!(mat, operator, basis)
    for (i,spinorI) in enumerate(basis)
        mat[:, i] = operator(spinorI)
    end
end
function getJtensor(L,β::Float64, theta; betaArray=zeros(Float64, L))
    if maximum(betaArray) == 0.0
         betas = [β^i for i in 0:L-2]
    else
        betas = betaArray
    end
    length(betas) != L-1 ? error("betaArray is wrong length.") : nothing

    jTensor = zeros(2^L, 2^L, Int(L*(L-1)/2))
    basis = getBasis(L)
    s2 = deepcopy(basis[1])
    s3 = deepcopy(basis[1])
    baseIndex = 0
    @views for k in 1:L-1
        for n in 1:L-k
            if betas[k] != 0.0
                operatorToMatrix!(jTensor[:,:,baseIndex + n], x->betas[k]*efficσiσj(n,n+k,x,s2,s3,theta), basis)
            end
        end
        baseIndex += L-k
    end
    return jTensor
end
function convertToArrays(tensor)
    smats = [sparse(tensor[:,:,i]) for i in 1:size(tensor)[3]]
    inds = [findnz(smat) for smat in smats]
    xs::Vector{Int64}, ys::Vector{Int64}, zs::Vector{Int64}, vs::Vector{Float64} = [], [], [], []
    for i in 1:size(tensor)[3]
        append!(xs, inds[i][1])
        append!(ys, inds[i][2])
        append!(zs, fill(i, length(inds[i][1])))
        append!(vs, inds[i][3])
    end
    return xs, ys, zs, vs
end
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
    #show(stdout,"text/plain", js[1:L])
    #show(stdout,"text/plain", jMat)
    return Hspace
end
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
@timeit to function getSpins!(ket, spinBasis, spinMatrix, ind)
   for i in 1:size(spinBasis)[1] #2^L
        for j in 1:size(spinBasis)[2] #L
            spinMatrix[j,ind] += abs2(ket[i])*spinBasis[i,j]
        end
    end
    return nothing
end


@timeit to function IsingefficU2(Hspace, hs, js, jArrays, hTensor; BCs="open") 
    efficientHam(Hspace, hs, js, jArrays, hTensor; BCs=BCs)
    for i in 1:size(Hspace)[1]
        Hspace[i,i] = exp(-im*Hspace[i,i])
    end
    return Hspace
    #return exp(-im.*efficientHam(Hspace, hs, js, jTensor, hTensor))
end

function efficU2(Hspace, hs, js, jTensor, hTensor) 
    return exp(-im.*Array(efficientHam(Hspace, hs, js, jTensor, hTensor)))
    #return exp(-im.*efficientHam(Hspace, hs, js, jTensor, hTensor))
end

density(mat::SparseMatrixCSC) = nnz(mat)/length(mat)
density(mat::Matrix) = length(findall(!iszero,mat))/length(mat)


@timeit to function newU1(L, ε)
    mat = zeros(2^L,2^L)
    basis = getBasis(L)
    answers = deepcopy(basis)    
    for i in 1:2^L
        for k in 1:L
            σx(k,basis[i],answers[k])
        end
        for k in 2:L
            answers[1] += answers[k] # This is where the allocations happen
        end 
        mat[:,i] .= answers[1]
    end
    if ε == 0.0
        return sparse(round.(exp(-im * mat .* (1-ε) * pi/2 ), digits=15))
    else
        return round.(exp(-im * mat .* (1-ε) * pi/2 ), digits=15)
    end
end

@timeit to function autocorrelator(spins, basis, eps, Ureal2, N)
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

        @timeit to "U2 mul" mul!(newKet,Ureal2,currentKet)
    
        getSpins!(newKet, basis, moreSpins, i)
        autoCor[:,i] = moreSpins[:,i] .* negOneSpins

        currentKet, newKet = newKet, currentKet
    end

    return autoCor, moreSpins
end
function gethsandjs(niters, L, J0, σj, σh)
    hs = zeros(niters, L)
    js = zeros(niters, Int(L*(L-1)/2))
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
@timeit to function effAvgAutoCor(niters, nperiods, spins, ε, J0, σj, σh; t=0.0, BCs="open")
    L = length(spins)
    Hspace = spzeros(ComplexF64, 2^L, 2^L)

    hs, js = gethsandjs(niters, L, J0, σj, σh)

    
    cors = zeros(L, nperiods+1, niters)
    finalCors=zeros(L, nperiods+1)
    allSpins = zeros(L, nperiods+1, niters)
    jIsingTensor = getIsingNNJtensor(L)
    #jTensor = getJtensor(L, 0.0, t)    # most of the time
    #jArrays = convertToArrays(jTensor) # also most of the tim
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

    return finalCors, allSpins[:,:,1]
end
@timeit to function avgLevelSpacings(niters, nperiods, spins, ε, J0, σj, σh; t=0.0, BCs="open")
    L = length(spins)
    Hspace = spzeros(ComplexF64, 2^L, 2^L)
    Hspace2 = Hermitian(spzeros(ComplexF64, 2^L, 2^L))

    hs, js = gethsandjs(niters, L, J0, σj, σh)
    rats = zeros(niters)
    jIsingTensor = getIsingNNJtensor(L)
    #jTensor = getJtensor(L, 0.0, t)    # most of the time
    #jArrays = convertToArrays(jTensor) # also most of the tim
    hTensor = getHtensor(L)            # more of the time

    basis = zeros(2^L, L)
    for i in 1:2^L
        basis[i,:] = Float64.(reverse(digits(i-1, base=2, pad=L)))
    end
    replace!(x->iszero(x) ? -1.0 : x, basis) #This basis is a matrix of the spins
    
    u1 = newU1(L, ε)
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

function lsrplotter(param, Lrange; BCs)

    if param == "eps"
        paramRange = range(0.0, 1.0, step=0.05)
    elseif param == "j"
        paramRange = 10 .^ collect(range(-2,1.5;step=0.2))
    elseif param == "sigH"
        paramRange = 10 .^ collect(range(-2,2;step=0.4))
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