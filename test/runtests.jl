using DTC, Test

@testset "src.jl tests" begin
    @testset "getIsingNNJtensor" begin
        @test getIsingNNJtensor(2) == [  1.0+0.0im  1.0+0.0im;
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
end

@testset "IO.jl tests" begin
    @testset "read/write vector" begin
        writeVector("testfile", [1.0, 2.0, 3.0])
        @test readVector("testfile") == [1.0, 2.0, 3.0]
        rm("testfile")
        @test_throws SystemError readVector("testfile")
    end
end