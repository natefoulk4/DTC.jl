using DTC, Test

@testset "basic stuff" begin
    @test getIsingNNJtensor(2) == [  1.0+0.0im  1.0+0.0im;
                                    -1.0+0.0im -1.0+0.0im;
                                    -1.0+0.0im -1.0+0.0im;
                                     1.0+0.0im  1.0+0.0im  ]
end