using DTC, Test

@testset "basic stuff" begin
    @test getIsingNNJtensor(2) == [  1.0+0.0im  1.0+0.0im;
                                    -1.0+0.0im -1.0+0.0im;
                                    -1.0+0.0im -1.0+0.0im;
                                     1.0+0.0im  1.0+0.0im  ]
    @test DTC.Rx(1, pi/2, [1.0+0.0im,  0.0+0.0im, 0.0+0.0im, 0.0+0.0im ], 
                            [1.0+0.0im,  0.0+0.0im, 0.0+0.0im, 0.0+0.0im ])\
                                ≈ [0.0+0.0im,  0.0+0.0im, 0.0+1.0im, 0.0+0.0im ]
end