using SymbolicControlSystems
using Test
using ControlSystems
using SymbolicControlSystems: s


@testset "SymbolicControlSystems.jl" begin
    @testset "To symbolic" begin
        @info "Testing To symbolic"


        pol = 2s^2+3s+4
        @test SymbolicControlSystems.expand_coeffs(pol) == [2,3,4]

        pol = 2s^3+3s^2+4*s
        @test SymbolicControlSystems.expand_coeffs(pol) == [2,3,4,0]

        @vars ω0
        sys = ControlSystems.DemoSystems.resonant(;ω0)
        ssys = Sym(sys)
        n, d = sp.fraction(simplify(ssys))
        @test SymbolicControlSystems.expand_coeffs(n) == [ω0^3]
        @test sum(SymbolicControlSystems.expand_coeffs(d) - [1, 0.5, ω0^2+0.25^2]) == 0



    end
end
