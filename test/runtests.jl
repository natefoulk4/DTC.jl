using DTC, Test, LinearAlgebra, Statistics

@testset "src.jl tests" begin
    @testset "getIsingJtensor" begin
        @test DTC.getIsingJtensor(2) == [  1.0+0.0im  1.0+0.0im;
                                        -1.0+0.0im -1.0+0.0im;
                                        -1.0+0.0im -1.0+0.0im;
                                        1.0+0.0im  1.0+0.0im  ]
    end
    @testset "Rx" begin
        @test DTC.Rx(1, pi/2, [1.0+0.0im,  0.0+0.0im, 0.0+0.0im, 0.0+0.0im ], 
                                [1.0+0.0im,  0.0+0.0im, 0.0+0.0im, 0.0+0.0im ]) ≈ 
                                [0.0+0.0im,  0.0+0.0im, 0.0+1.0im, 0.0+0.0im ]
    end
    @testset "getHtensor" begin
        res = DTC.getHtensor(2)
        @test res == [-1.0 -1.0; -1.0 1.0; 1.0 -1.0; 1.0 1.0]
        l = rand(1:5)
        res = DTC.getHtensor(l)
        @test size(res) == (2^l, l)
    end
    @testset "levelspacing" begin
        @test levelspacing([2.0, 4.0, 3.0, 1.0]) == levelspacing([1.0, 2.0, 3.0, 4.0])
        @test levelspacing([2.0, 4.0, 3.0, 1.0]) == 1.0
        @test levelspacing([1.0, 3.0, 4.0, 5.0]) == 0.75
        rng = rand(Float64,4)
        @test levelspacing(rng) ≈ levelspacing(rng .+ 0.5)
        @test_throws MethodError levelspacing([1.0 3.0 4.0 5.0])
        @test_throws MethodError levelspacing([1.0+0.0im, 2.0+0.0im, 3.0+0.0im, 4.0+0.0im])
    end
    @testset "getBasis" begin 
        @test DTC.getBasis(2) == complex.([-1 -1; -1 +1; +1 -1; +1 +1])
        @test DTC.getBasis(3) == complex.([-1 -1 -1; -1 -1 +1; -1 +1 -1; -1 +1 +1;
                                           +1 -1 -1; +1 -1 +1; +1 +1 -1; +1 +1 +1])
    end
    @testset "getJsAndHs" begin
        testParams = DTC.getJsAndHs(10, 5, 1.0, 0.0, "discrete", 3.0, 0.0, "normal")
        @test size(testParams[1]) == (10, 5) == size(testParams[2])
        @test all(map(isequal(1.0), testParams[1]))
        @test all(map(isequal(3.0), testParams[2]))
    end
    @testset "getKet" begin
        @test DTC.getKet([1,0,1]) == complex.([0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0] )
        @test DTC.getKet([0,0,1]) == complex.([0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) 
        @test DTC.getKet([1,1,1]) == complex.([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0])
        @test_throws MethodError DTC.getKet([1.0, 0.0, 1.0])
    end
    @testset "binary to base10" begin
        @test binToBase10([1,0,1]) == 5
        @test binToBase10([0]) == 0
        @test_throws MethodError binToBase10(10)
    end
    @testset "base10 to binary" begin
        @test base10ToBin(8,3) == [1,1,1]
        @test base10ToBin(1,4) == [-1,-1,-1,-1]
    end
    @testset "efficientHam" begin
        @test DTC.efficientHam(zeros(ComplexF64, (4,4)), [-1.0, 1.0], [3.0], [1.0 -1.0 -1.0 1.0]', [-1.0 -1.0 1.0 1.0; -1.0 1.0 -1.0 1.0]'; theta=0.0, BCs="open") == convert.(ComplexF64, [3.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0; 0.0 0.0 -5.0 0.0; 0.0 0.0 0.0 3.0])
    end
    @testset "getSpins!" begin
        spins = [0.0; 0.0]
        DTC.getSpins!([0.0, 1.0, 0.0, 0.0], DTC.getBasis(2), spins, 1)
        @test spins[:,1] == [-1.0, 1.0]
        spins = [0.0; 0.0; 0.0]
        DTC.getSpins!(normalize([0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0]), DTC.getBasis(3), spins, 1)
        @test spins[:,1] ≈ [1.0, 0.0, 1.0]
    end
    @testset "autocorrelator" begin
        @test autocorrelator([1,0,1,0],DTC.getBasis(4),0.0, DTC.IsingefficU2(zeros(ComplexF64, (16,16)), [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0], DTC.getIsingJtensor(4), DTC.getHtensor(4); t=0.0), 100)[:,end] ≈ DTC.getKet([1,0,1,0])
        @test autocorrelator([1,0,1,0],DTC.getBasis(4),0.0, DTC.IsingefficU2(zeros(ComplexF64, (16,16)), [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0], DTC.getIsingJtensor(4), DTC.getHtensor(4); t=0.0), 100)[:,end-1] ≈ DTC.getKet([0,1,0,1])

        @test autocorrelator([1,0,0,1],DTC.getBasis(4),0.0, DTC.IsingefficU2(zeros(ComplexF64, (16,16)), [-1.0, 1.0, 1.0, -1.0], [1.0, 1.0, 1.0, 1.0], DTC.getIsingJtensor(4), DTC.getHtensor(4); t=0.0), 100)[:,end] ≈ DTC.getKet([1,0,0,1])

        finalKet = autocorrelator([1,0,1,0,0,0,0],DTC.getBasis(7), 0.1, DTC.IsingefficU2(zeros(ComplexF64, (2^7,2^7)), [-0.3, 0.4, -0.5, 0.6, -0.7, 0.8, -0.9], [1.1, 1.05, -1.02, 1.01, -0.95, -0.9, -0.98], DTC.getIsingJtensor(7), DTC.getHtensor(7); t=0.0), 100)[:,end] 
        spins = zeros(7,1)
        DTC.getSpins!(finalKet, DTC.getBasis(7), spins, 1)
        spins[1] ≈ 0.9586187956175708
    end
    @testset "effAvgAutoCor" begin
        @test isapprox(mean(effAvgAutoCor(2000, 300, [1,0,1,0]; ε=0.10, J0=1.0, σj=0.10, J_dist="discrete", H_dist="normal", σh=2.0, t=0.0)[1][:,end]),   0.8177; atol=0.02)
        @test isapprox(mean(effAvgAutoCor(2000, 300, [1,0,1,0]; ε=0.10, J0=1.0, σj=0.10, J_dist="discrete", H_dist="normal", σh=2.0, t=0.0, diagonalization=true)[1][:,end]),   0.8177; atol=0.02)
        @test effAvgAutoCor(1,200, [1,0,1,0]; ε=0.1, J0=pi/4, σj=0.0, h0=0.05, σh=0.0, t=1.0, num_H2I=64, J_dist="normal", H_dist="normal", BCs="open")[1][1,end] ≈ 0.9406572121389324
        @test effAvgAutoCor(1,200, [1,0,1,0]; ε=0.1, J0=pi/4, σj=0.0, h0=0.05, σh=0.0, t=1.0, num_H2I=64, J_dist="normal", H_dist="normal", BCs="open", diagonalization=true)[1][1,end] ≈ 0.9406572121389324
    end
    @testset "U1" begin
        @test U1(2,0.0) == complex.([0.0 0.0 0.0 -1.0; 
                              0.0 0.0 -1.0 0.0;
                              0.0 -1.0 0.0 0.0;
                              -1.0 0.0 0.0 0.0])
        @test matrix_density(U1(3, 0.0)) == 0.125
        @test matrix_density(U1(3, 0.05)) == 1.0
        @test matrix_density(U1(4, 0.0)) == 0.0625
        @test matrix_density(U1(4, 0.05)) == 1.0
    end
end

@testset "IO.jl tests" begin
    @testset "read/write vector" begin
        writeArray("testfile", [1.0, 2.0, 3.0])
        @test readArray("testfile") == [1.0, 2.0, 3.0]
        rm("testfile")
        @test_throws SystemError readArray("testfile")
    end
end