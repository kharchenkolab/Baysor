using Baysor, Test

@testset "cli" begin
    @testset "Help" begin
        redirect_stdout(open("/dev/null", "w")) do; Baysor.run_cli(["--help"]); end
        @test 1 == 1
    end
end
