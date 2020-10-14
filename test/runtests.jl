using SubglacialPlumes
using Test


#params = Params( )
@testset "ExampleTests" begin
    #2x + 3y
    @test f(2,1) == 7
    @test f(2,3) == 13
end
