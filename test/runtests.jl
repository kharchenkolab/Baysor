using Baysor, Test

@testset "cli" begin
    @testset "Help" begin
        Baysor.run_cli(["--help"])
        @test 1 == 1
    end
end
