using DTC, Test, LinearAlgebra, Statistics

@testset "src.jl tests" begin
    @testset "getIsingNNJtensor" begin
        @test DTC.getIsingNNJtensor(2) == [  1.0+0.0im  1.0+0.0im;
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
        @test_throws MethodError levelspacing([1.0 3.0 4.0 5.0])
        @test_throws MethodError levelspacing([1.0+0.0im, 2.0+0.0im, 3.0+0.0im, 4.0+0.0im])
    end
    @testset "getBasis" begin 
        @test DTC.getBasis(2) == complex.([-1 -1; -1 +1; +1 -1; +1 +1])
        @test DTC.getBasis(3) == complex.([-1 -1 -1; -1 -1 +1; -1 +1 -1; -1 +1 +1;
                                           +1 -1 -1; +1 -1 +1; +1 +1 -1; +1 +1 +1])
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
        @test DTC.efficientHam(zeros(ComplexF64, (4,4)), [-1.0, 1.0], [3.0], [1.0 -1.0 -1.0 1.0]', [-1.0 -1.0 1.0 1.0; -1.0 1.0 -1.0 1.0]'; BCs="open") == convert.(ComplexF64, [3.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0; 0.0 0.0 -5.0 0.0; 0.0 0.0 0.0 3.0])
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
        @test minimum(autocorrelator([1,0,1,0],DTC.getBasis(4),0.0, DTC.IsingefficU2(zeros(ComplexF64, (16,16)), [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0], DTC.getIsingNNJtensor(4), DTC.getHtensor(4)), 100)[1][1,:][1:2:end]) ≈ 1.0
        @test minimum(autocorrelator([1,0,0,1],DTC.getBasis(4),0.0, DTC.IsingefficU2(zeros(ComplexF64, (16,16)), [-1.0, 1.0, 1.0, -1.0], [1.0, 1.0, 1.0, 1.0], DTC.getIsingNNJtensor(4), DTC.getHtensor(4)), 100)[1][1,:][1:2:end]) ≈ 1.0
        @test autocorrelator([1,0,1,0,0,0,0],DTC.getBasis(7), 0.1, DTC.IsingefficU2(zeros(ComplexF64, (2^7,2^7)), [-0.3, 0.4, -0.5, 0.6, -0.7, 0.8, -0.9], [1.1, 1.05, -1.02, 1.01, -0.95, -0.9, -0.98], DTC.getIsingNNJtensor(7), DTC.getHtensor(7)), 100)[1][1,:][end] ≈ 0.9586187956175708
    end
    @testset "effAvgAutoCor" begin
        isapprox(mean(effAvgAutoCor(2000, 300, [1,0,1,0], 0.10, 1.0, 0.10, 2.0)[1][:,end]),   0.8177; atol=0.02)
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
        writeVector("testfile", [1.0, 2.0, 3.0])
        @test readVector("testfile") == [1.0, 2.0, 3.0]
        rm("testfile")
        @test_throws SystemError readVector("testfile")
    end
end