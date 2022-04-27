
using LinearAlgebra, Statistics, Distributions, TimerOutputs


"
binToBase10(x::Vector{Int})
Convert a binary vector into a base-10 number."
function binToBase10(x)
    replace!(y->isequal(-1.0, y) ? 0.0 : y, x)
    L = length(x)
    basis = [2^(L-i) for i in 1:L]
    return (basis⋅x)
end

"
base10ToBin(x)
Convert a a base-10 number into a binary vector of length L"
function base10ToBin(x, L)
    return replace!(x->iszero(x) ? -1.0 : x, complex.(reverse(digits(x-1, base=2, pad=L))))
end

"
    Rx(n, φ, spinKet, newKet)
Operate on ``spinKet`` (length ``2^L``) with the rotation unitary R^x_n(φ) (on the ``n``th qubit). Return ``newKet``"
function Rx(n::Int, φ::Real, spinKet::Vector{ComplexF64}, newKet::Vector{ComplexF64})
    L = Int(log2(length(spinKet)))
    stride =  2^(L-n)   # 1 if n=4, 2 if n=3, 4 if n=2, 8 if n = 1
    stride2 = 2^(n-1)   # 8 if n=4, 4 if n=3, 2 if n=2, 1 if n = 1

    cos1 = cos(φ)
    sin1 = sin(φ)

    for i in 1:stride2
        for j in 1:stride
            #println(2*(i-1)*stride+j, 2*(i-1)*stride+j + stride)
            newKet[2*(i-1)*stride+j] = cos1*spinKet[2*(i-1)*stride+j] + 
                                            im*sin1*spinKet[2*(i-1)*stride + j + stride]
            newKet[2*(i-1)*stride+j + stride]  = im*sin1*spinKet[2*(i-1)*stride+j] +
                                            cos1*spinKet[2*(i-1)*stride+j + stride]
        end
    end

    return newKet
end

"
    Rz(n, φ, spinKet, newKet)
Operate on ``spinKet`` (length ``2^L``) with the rotation unitary R^z_n(φ) (on the ``n``th qubit). Return ``newKet``"
function Rz(n::Int, φ::Real, spinKet::Vector{ComplexF64}, newKet::Vector{ComplexF64})
    L = Int(log2(length(spinKet)))
    stride =  2^(L-n)   # 1 if n=4, 2 if n=3, 4 if n=2, 8 if n = 1
    stride2 = 2^(n-1)   # 8 if n=4, 4 if n=3, 2 if n=2, 1 if n = 1
    exp1 = exp(im*φ)
    exp2 = exp(-im*φ)

    for i in 1:stride2
        for j in 1:stride
            newKet[2*(i-1)*stride+j] = exp1*spinKet[2*(i-1)*stride+j]
            newKet[2*(i-1)*stride+j + stride]  = exp2*spinKet[2*(i-1)*stride+j + stride]
        end
    end

    return newKet
end

"
    getIsingJtensor(L)
Return ``2^L x L`` matrix, where the ``n``th column contains the approriate signs for the diagonal of the J\\_n matrix. _Ising model only._"
@timeit to function getIsingJtensor(L)
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
    getHeisenCoords(L)
    Return the coordinates of the heisenberg non-diagonal elements of the hamiltonian."
function getHeisenCoords(L)
    coords = zeros(Int64, (2, 2^(L-2), L))

    for n in 1:L-1
        k=1
        for j in 1:2^(L-n-1)
            for i in 1:2^(n-1) 
                stride = 2^(n+1)
                coords[1,k,n] = 2^(n-1) + (j-1)*stride + i
                coords[2,k,n] = 2^(n) + (j-1)*stride + i
                k += 1
            end
        end
    end

    for i in 1:2^(L-2)
        coords[1,i,L] = 2*i
        coords[2,i,L] = 2*i + 2^(L-1) - 1
    end

    return coords
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
@timeit to function levelspacing(vals::Vector{<:Real})
    sortedvals = sort(vals)
    for i in 1:length(sortedvals)-1
        sortedvals[i] = sortedvals[i+1] - sortedvals[i]
    end
    for i in 1:length(sortedvals)-2
        if sortedvals[i] > sortedvals[i+1]
            sortedvals[i] = sortedvals[i+1]/sortedvals[i]
        elseif sortedvals[i+1] > sortedvals[i] 
            sortedvals[i] = sortedvals[i]/sortedvals[i+1]
        else
            sortedvals[i] = 1.0
        end
    end
    return mean(sortedvals[1:end-2])
end

"
    efficientHam(Hspace, hs, js, jMat, hTensor; tCoords=nothing, theta, BCs)
Multiply each J and H columns of ``jMat`` and ``hTensor`` by the appropriate prefactors contained in ``js`` and ``hs``. Add them together and assign them to the diagonal of ``Hspace``, overwriting ``Hspace`` completely. Return ``Hspace``. keyword: ``BCs={'open', 'periodic'}``"
@timeit to function efficientHam(Hspace, hs, js, jMat, hTensor; tCoords=nothing, theta, BCs)
    L = Int(log2(size(Hspace)[1]))
    Hspace .= 0.0 + 0.0im
    if BCs == "open"
        for n in 1:L-1 
            for j in 1:2^L
                Hspace[j,j] += (js[n] * jMat[j,n])
            end
            if theta != 0.0
                for k in 1:2^(L-2)
                    Hspace[tCoords[1,k,n], tCoords[2,k,n]] += js[n]*2.0*theta
                    Hspace[tCoords[2,k,n], tCoords[1,k,n]] += js[n]*2.0*theta
                end
            end 
        end
    elseif BCs == "periodic"
        for n in 1:L
            for j in 1:2^L
                Hspace[j,j] += (js[n] * jMat[j,n])
            end
            if theta != 0.0
                for k in 1:2^(L-2)
                    Hspace[tCoords[1,k,n], tCoords[2,k,n]] += js[n]*2.0*theta
                    Hspace[tCoords[2,k,n], tCoords[1,k,n]] += js[n]*2.0*theta
                end
            end 
        end
    end

    for j in 1:L, n in 1:2^L
        Hspace[n,n] += hs[j] * hTensor[n,j]
    end
    #show(stdout, "text/plain", round.(Float64.(Array(Hspace)),digits=2))
    return Hspace
end

"
    getKet(spinArray)
Given a (length ``L``) array of spins (must be a pure state), find the appropriate base ket (length ``2^L``) of Hilbert space."
function getKet(spinArray::Vector{Int})
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
@timeit to function IsingefficU2(Hspace, hs, js, jMat, hTensor; tCoords=nothing, t=0.0, BCs="open", num_H2I=0) 
    efficientHam(Hspace, hs, js, jMat, hTensor; tCoords=tCoords, theta=t, BCs=BCs)
    if t == 0.0
        for i in 1:size(Hspace)[1]
            Hspace[i,i] = exp(-im*Hspace[i,i])
        end
        return Hspace
    elseif num_H2I % 2 == 0
        denom = num_H2I > 0 ? num_H2I : 1
        return sparse(exp(Array(-im.*Hspace ./ denom)))
    else
        error("theta and num_H2I are incorrectly specified!!")
    end
end

"
    matrix_density(mat)
Given a matrix ``mat``, return the density of the matrix, which is the ratio of nonzero entries to the total number of matrix entries."
matrix_density(mat::Matrix) = length(findall(!iszero,mat))/length(mat)

"
    autocorrelator(spins, basis, eps, U2, N)
Simulate the dynamics of a DTC for a given instance of disorder. Take initial state of ``spins`` (length ``L``) and ``basis`` (size ``2^L x L``). Apply Floquet unitaries (using ``eps``, the perturbation of a π pulse, and ``U2``) ``N`` times. Return a matrix of all the kets (``2^L`` x ``N+1``)."
@timeit to function autocorrelator(spins, basis, eps, U2, N; num_H2I=0, d=false)
    initKet = getKet(spins)
    L = length(spins)

    @timeit to "initializing kets" begin
        allKets = zeros(ComplexF64, 2^L, N+1)
        currentKet = deepcopy(initKet)
        interKet = zeros(ComplexF64, 2^L)
        newKet = zeros(ComplexF64, 2^L)
    end


    @timeit to "Rx,Rz, and U2" begin
        if d 
            vals = zeros(ComplexF64, 2^L)
            bigU = zeros(ComplexF64, (2^L,2^L))
            u1 = U1(L, eps)
            if num_H2I > 0
                bigU = u1
                for k in 1:Int(num_H2I/2)
                    bigU = U2*bigU
                    bigU = Uz(L, -pi/2)*bigU
                    bigU = U2*bigU
                    bigU = Uz(L, pi/2)*bigU
                end
            else
                bigU = U2*u1
            end

            @timeit to "eigen" factorization = eigen(bigU)
            vals .= factorization.values
            currentKet = transpose(conj.(factorization.vectors)) * currentKet
            allKets[:,1] .= currentKet
            for i in 2:N+1
                currentKet .=  vals .* currentKet
                allKets[:,i] .= currentKet
            end

            @timeit to "matmul at end" allKets = factorization.vectors * allKets
        else 
            allKets[:,1] = currentKet
            for i in 2:N+1
                for k in 1:L
                    #mul!(interKet, U1s[k], currentKet)
                    Rx(k, pi/2*(1-eps), currentKet, interKet)
                    interKet, currentKet = currentKet, interKet
                end
                if num_H2I > 0
                    for k in  1:Int(num_H2I/2)
                        mul!(interKet,U2,currentKet)
                        for j in 1:2:L
                            Rz(j, -pi/2, interKet, currentKet)
                            interKet, currentKet = currentKet, interKet
                        end
                        mul!(currentKet,U2,interKet)
                        for j in 1:2:L
                            Rz(j, pi/2, currentKet, interKet)
                            interKet, currentKet = currentKet, interKet
                        end
                    end
                    newKet .= currentKet
                else
                    mul!(newKet,U2,currentKet)
                end
                allKets[:,i] = newKet
                currentKet, newKet = newKet, currentKet
            end
        end
    end

    allKets = abs2.(allKets)

    return allKets
end

"
    getJsAndHs(niters, L, J0, σJ, J_dist, h0, σh, h_dist)
Sample distribution ``niter*L`` times. Return Hmatrix (``niters x L``) and J matrix (``niters x L``)"
function getJsAndHs(niters, L, J0, σj, J_dist, h0, σh, h_dist)
    hs = zeros(niters, L)
    js = zeros(niters, L)
    M = 8
    discrete_Jrange = J0 .+ (σj .* cos.(pi*collect(0:M-1)/(M-1)))

    for i in eachindex(hs)
        if σh > 0.0
            if lowercase(h_dist) == "normal"
                hs[i] = rand(Normal(h0, σh))
            elseif lowercase(h_dist) == "uniform"
                hs[i] = rand(Uniform(h0-σh, h0+σh))
            else
                error("No supported distribution for onsite disorder given.")
                return nothing
            end
        else
            hs[i] = h0
        end
    end

    for i in eachindex(js)
        if σj > 0.0
            if lowercase(J_dist) == "normal"
                js[i] = rand(Normal(J0, σj))
            elseif lowercase(J_dist) == "uniform"
                js[i] = rand(Uniform(J0-σj, J0+σj))
            elseif lowercase(J_dist) == "discrete"
                js[i] = rand(discrete_Jrange)
            else
                error("No supported distribution for nearest neighbor disorder given.")
                return nothing
            end
        else
            js[i] = J0
        end
    end

    return js, hs
end

@timeit to function U1(L, ε)
    mat1 = zeros(ComplexF64, (2^L, 2^L))
    mat2 = UniformScaling(1)
    ket1 = zeros(ComplexF64, 2^L)
    ket2 = zeros(ComplexF64, 2^L)
    for n in 1:L
        for i in 1:2^L
            ket1 .= 0.0+0.0im
            ket1[i] = 1.0+0.0im
            Rx(n, (1-ε)*pi/2, ket1, ket2)
            mat1[:,i] = ket2
        end
        mat2 = mat1 * mat2 
    end
    return round.(mat2, digits=15)
end

@timeit to function Uz(L, φ)
    mat1 = zeros(ComplexF64, (2^L, 2^L))
    mat2 = UniformScaling(1)
    ket1 = zeros(ComplexF64, 2^L)
    ket2 = zeros(ComplexF64, 2^L)
    for n in 1:2:L
        for i in 1:2^L
            ket1 .= 0.0+0.0im
            ket1[i] = 1.0+0.0im
            Rz(n, φ, ket1, ket2)
            mat1[:,i] = ket2
        end
        mat2 = mat1 * mat2 
    end
    return round.(mat2, digits=15)
end

"
    effAvgAutoCor(niters, nperiods, spins; ε, J0, σJ, J_dist, h0=0.0, σH, H_dist, t=0.0, BCs='open', verbose=false)
Exact same as autocorrelator, except average over ``niters`` simulations. Return ``cors`` and ``allSpins``, which are both ``2^L x L`` matrices "
@timeit to function effAvgAutoCor(niters, nperiods, spins; ε, J0, σj, J_dist="normal", h0=0.0, σh, H_dist="normal", t=0.0, BCs="open", num_H2I=0, diagonalization=true, verbose=false)
    if BCs != "open" && BCs != "periodic"
        error("Boundary conditions must be 'open' or 'periodic' !")
    end
    #nothing
    if t == 0.0
        numH2I = 0
    end
    L = length(spins)
    Hspace = spzeros(ComplexF64, 2^L, 2^L)
    negOneSpins = replace(spins, 0 => -1)

    js, hs =getJsAndHs(niters, L, J0, σj, J_dist, h0, σh, H_dist)

    
    cors = zeros(L, nperiods+1)
    allSpins = zeros(L, nperiods+1)
    @timeit to "initializing allKets" allKets = zeros(2^L, nperiods+1)
    jIsingTensor = getIsingJtensor(L)
    thetaCoords = t == 0.0 ? nothing : getHeisenCoords(L)
    hTensor = getHtensor(L)            # more of the time
    basis = getBasis(L)

    for i in 1:niters
        u2 = IsingefficU2(Hspace,  hs[i,:] ,  js[i,:], jIsingTensor, hTensor; tCoords=thetaCoords, t=t, BCs=BCs, num_H2I=num_H2I)
        allKets += autocorrelator(spins, basis, ε, u2, nperiods; num_H2I=num_H2I, d=diagonalization)
        if i % (niters/10) == 0  && verbose == true
            println("Finished ",i,"th iteration")
        end
    end

    @timeit to "averaging over iters" allKets = allKets ./ niters
    @timeit to "getting spins" allSpins = basis' * allKets
    @timeit to "getting cors" cors = allSpins .* negOneSpins

    return cors, allSpins
end

"   
    getBasis(L)
Return a ``2^L x L`` matrix where the ``i``th row corresponds to the ``i``th spin ket (+1 and -1)."
function getBasis(L)
    basis = zeros(2^L, L)
    for i in 1:2^L
        basis[i,:] = base10ToBin(i, L)
    end
    return basis
end

"
    avgLevelSpacings(niters, spins; ε, J0, σJ, J_dist, h0=0.0, σH, H_dist, t=0.0, BCs='open')
Calculate ``niters`` of the level spacing ratios. Return the average LSR overall."
@timeit to function avgLevelSpacings(niters, spins; ε, J0, σj, J_dist, h0=0.0, σh, H_dist, t=0.0, BCs="open")
    L = length(spins)
    Hspace = zeros(ComplexF64, 2^L, 2^L)
    Hspace2 = Hermitian(zeros(ComplexF64, 2^L, 2^L))

    js, hs =getJsAndHs(niters, L, J0, σj, J_dist, h0, σh, H_dist)
    rats = zeros(niters)
    jIsingTensor = getIsingJtensor(L)
    hTensor = getHtensor(L)            # more of the time

    u1 = U1(L, ε)
    noconverge=0

    for i in 1:niters
        @timeit to "matmul" Hspace2 = Array(u1*IsingefficU2(Hspace,hs[i,:],js[i,:],jIsingTensor,hTensor, BCs=BCs))
        @timeit to "eigs" vals = eigvals!(Hspace2)
        
        rats[i] = levelspacing( mod.(Real.(round.(log.(Complex.(vals)) .* im,digits=8)), 2*pi) )
        if i % (niters/10) == 0
            #println("Finished $i th iteration out of $niters")
        end
    end

    return mean(rats)
end

"
    LsrsOverParamRange(npoints, niters, param, Lrange; BCs)
Calculate LSRs for ``npoints`` points over a parameter range (``param={'eps','j','sigH'}``), averaging each point over ``niters`` times. Done for several values of L (``Lrange``). Boundary conditions (``BCs={'open','periodic'}``) optional."
function LsrsOverParamRange(npoints, niters, param, Lrange; BCs)

    if param == "eps"
        paramRange = range(0.0, 1.0, length=npoints)
    elseif param == "j"
        paramRange = 10 .^ collect(range(-2,1.5; length=npoints))
    elseif param == "sigH"
        paramRange = 10 .^ collect(range(-2,2; length=npoints))
    else
        error("Parameter not correctly specified!")
    end
    n = length(paramRange)
    lsrs = zeros(length(Lrange), n)

    niter = niters
    J0 = pi/4
    σJ = pi/8
    J_dist="discrete"
    h0 = 0.0
    σH  = pi/50
    H_dist="normal"
    ε = 0.1

    i = 1

    for l in Lrange
        println("Starting L=",l)
        for (j,paramValue) in enumerate(paramRange) 
            if param == "eps"
                ε = paramValue
            elseif param == "j"
                J0 = paramValue
            elseif param == "sigH"
                σH = paramValue
            end
            init = rand([0,1], l)
            lsrs[i,j] = avgLevelSpacings(niter, init; ε=ε, J0=J0, σj=σJ, J_dist=J_dist, h0=h0, σh=σH, H_dist=H_dist, BCs=BCs)
            println("Finished $j th data point out of $n ($niter iterations each point)")
        end
        i += 1
    end

    return paramRange, lsrs'
end