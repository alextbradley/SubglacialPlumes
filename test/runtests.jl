using SubglacialPlumes
using Test


#params = Params( )
@testset "ExampleTests" begin
    #2x + 3y
    @test SubglacialPlumes.f(2,1) == 7
    @test SubglacialPlumes.f(2,3) == 13
end
