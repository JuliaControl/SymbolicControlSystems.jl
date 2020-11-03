using SymbolicControlSystems
using Test
using ControlSystems
using SymbolicControlSystems: s, z


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

        @testset "Tustin and C-code" begin
            @info "Testing Tustin and C-code"
            
            @vars J c
            G = tf(1.,[J^2,c,1])
            Gn = tf(1.,[1,1,1])
            Gt = tustin(G, 0.1)
            @test_throws ErrorException sym2num(Gt, J=>1, c=>1)
            Gtn = sym2num(Gt, 0.1, J=>1, c=>1)
            @test Gtn ≈ c2d(Gn, 0.1, :tustin)
            code = SymbolicControlSystems.ccode(Gt)
            @test occursin("static double u[3] = {0};", code)
            @test occursin("static double y[3] = {0};", code)
            @test occursin("double ui, double J, double c", code)
            @test occursin("""#include <stdio.h>

double transfer_function(double ui, double J, double c) {
    static double u[3] = {0};
    static double y[3] = {0};
    int i;
    for (i=2; i > 0; --i) {
        u[i] = u[i-1];
    }
    u[0] = ui;
    for (i=2; i > 0; --i) {
        y[i] = y[i-1];
    }
    y[0] = 0;
    y[0] += (1.0/(400.0*pow(J, 2) + 20.0*c + 1.0))*u[0];
    y[0] += (2.0/(400.0*pow(J, 2) + 20.0*c + 1.0))*u[1];
    y[0] += (1.0/(400.0*pow(J, 2) + 20.0*c + 1.0))*u[2];
    y[0] += ((800.0*pow(J, 2) - 2.0)/(400.0*pow(J, 2) + 20.0*c + 1.0))*y[1];
    y[0] += ((-400.0*pow(J, 2) + 20.0*c - 1.0)/(400.0*pow(J, 2) + 20.0*c + 1.0))*y[2];
    return y[0];
}""", code)
        end
            

    end
end
