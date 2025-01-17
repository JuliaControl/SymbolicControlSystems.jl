using SymbolicControlSystems
using Test
using ControlSystemsBase
s = SymbolicControlSystems.s
z = SymbolicControlSystems.z
import Symbolics
import Symbolics: Num
using LinearAlgebra
Base.active_repl.options.hint_tab_completes = false # This messes with sympy https://discourse.julialang.org/t/sympy-makes-repl-to-stuck/124814/6


macro test_both(G,sym)
    quote
        if isdiscrete($(esc(G))) 
            @test Sym($(esc(G))) == $(esc(sym))
            @test $(esc(G)) == tf($(esc(sym)), $(esc(G)).Ts)
        else
            @test Sym($(esc(G))) == $(esc(sym))
            @test $(esc(G)) == tf($(esc(sym)))
            @test tf(ss($(esc(G)))) == tf($(esc(sym)))
        end
    end
end

macro test_both_symb(G,sym)
    quote
        if isdiscrete($(esc(G))) 
            @test isequal(to_num($(esc(G))), $(esc(sym)))
            @test isequal($(esc(G)), tf($(esc(sym)), $(esc(G)).Ts))
        else
            @test isequal(to_num($(esc(G))), $(esc(sym)))
            @test $(esc(G)) == tf($(esc(sym)))
            # @test isequal(tf(ss($(esc(G)))), tf($(esc(sym)))) # ss -> tf calls eigvals which fails on Num
        end
    end
end


@testset "SymbolicControlSystems.jl" begin
    @testset "To symbolic" begin
        @info "Testing To symbolic"


        pol = 2s^2+3s+4
        @test SymbolicControlSystems.expand_coeffs(pol, s) == [2,3,4]

        pol = 2s^3+3s^2+4*s
        @test SymbolicControlSystems.expand_coeffs(pol, s) == [2,3,4,0]

        pol = 2s^3+4*s
        @test SymbolicControlSystems.expand_coeffs(pol, s) == [2,0,4,0]

        pol = 2z^2+3z+4
        @test SymbolicControlSystems.expand_coeffs(pol, z) == [2,3,4]

        pol = 2z^3+3z^2+4*z
        @test SymbolicControlSystems.expand_coeffs(pol, z) == [2,3,4,0]

        pol = 2z^3+4*z
        @test SymbolicControlSystems.expand_coeffs(pol, z) == [2,0,4,0]

        @syms ω0
        sys = ControlSystemsBase.DemoSystems.resonant(;ω0)
        ssys = Sym(sys)
        n, d = sp.fraction(simplify(ssys))
        @test SymbolicControlSystems.expand_coeffs(n, s) == [1.0*ω0^3]
        @test sum(SymbolicControlSystems.expand_coeffs(d, s) - [1, 0.5, ω0^2+0.25^2]) == 0
    end

    @testset "SymPy: Sym -> tf and vice versa" begin
        @info "SymPy: Testing Sym -> tf and vice versa"

        @syms a b c
        @test_both tf([a], [b, c])  a/(b*s + c)
        @test_both tf([a], [b, c], 0.1)  a/(b*z + c)

        @test_both tf([a], [b, c, 1])  a/(b*s^2 + c*s + 1)
        @test_both tf([a], [b, c, 1], 0.1)  a/(b*z^2 + c*z + 1)

        # Test that denominator is still monic
        @test string(Sym(tf([a], [1, c, b]))) == string(a/(s^2 + c*s + b))
        @test string(Sym(tf([a], [1, c, b], 0.1))) == string(a/(z^2 + c*z + b))


        # Test that denominator is still monic
        @test string(Sym(tf([b, a], [1, c, 1]))) == string((b*s + a)/(s^2 + c*s + 1))
        @test string(Sym(tf([b, a], [1, c, 1], 0.1))) == string((b*z + a)/(z^2 + c*z + 1))


        # Test that denominator is still monic
        @test_both tf([a], [1, c, b])  a/(s^2 + c*s + b)
        @test_both tf([a], [1, c, b], 0.1)  a/(z^2 + c*z + b)


        # Test that denominator is still monic
        @test_both tf([b, a], [1, c, 1])  (b*s + a)/(s^2 + c*s + 1)
        @test_both tf([b, a], [1, c, 1], 0.1)  (b*z + a)/(z^2 + c*z + 1)
        
    end


    @testset "sym2num" begin
        @info "Testing sym2num"
        @syms a b c
        G = tf([a], [b, c]) 
        @test sym2num(G, a=>1, b=>1, c=>1) == tf(1,[1,1])
        @test sym2num(G, a=>1, b=>2, c=>1) == tf(1,[2,1])
        @test sym2num(G, a=>1, b=>1, c=>2) == tf(1,[1,2])

        G = tf([b, a], [1, c, 1])
        @test sym2num(G, a=>1, b=>1, c=>1) == tf([1,1],[1,1,1])
        @test sym2num(G, a=>1, b=>2, c=>1) == tf([2,1],[1,1,1])
        @test sym2num(G, a=>1, b=>1, c=>2) == tf([1,1],[1,2,1])

        G = tf([a], [b, c], 0.1) 
        @test_throws ErrorException sym2num(G, a=>1, b=>1, c=>1) == tf(1,[1,1], 0.1)
        
        @test sym2num(G, 0.1, a=>1, b=>1, c=>1) == tf(1,[1,1], 0.1)
        @test sym2num(G, 0.1, a=>1, b=>2, c=>1) == tf(1,[2,1], 0.1)
        @test sym2num(G, 0.1, a=>1, b=>1, c=>2) == tf(1,[1,2], 0.1)

        G = tf([b, a], [1, c, 1], 0.1)
        @test sym2num(G, 0.1, a=>1, b=>1, c=>1) == tf([1,1],[1,1,1], 0.1)
        @test sym2num(G, 0.1, a=>1, b=>2, c=>1) == tf([2,1],[1,1,1], 0.1)
        @test sym2num(G, 0.1, a=>1, b=>1, c=>2) == tf([1,1],[1,2,1], 0.1)
        
    end

    @testset "Tustin and C-code" begin
        @info "Testing Tustin and C-code"
        
        @syms J c
        G = tf(1.,[J^2,c,1])
        Gs = Sym(G)
        Gn = tf(1.,[1,1,1])
        Gt = tustin(G, 0.1)
        @test denvec(Gt)[1][1] == 1 # test that it's monic
        @test_throws ErrorException sym2num(Gt, J=>1, c=>1)
        Gtn = sym2num(Gt, 0.1, J=>1, c=>1)
        @test minreal(Gtn) ≈ c2d(Gn, 0.1, :tustin)
        code = SymbolicControlSystems.ccode(Gt)
        @test occursin("static double u[3] = {0};", code)
        @test occursin("static double y[3] = {0};", code)
        @test occursin("double ui, double J, double c", code)
        
        code = SymbolicControlSystems.ccode(ss(Gt), cse=false)
        @test occursin("double u, double J, double c", code)
       
        code = SymbolicControlSystems.ccode(ss(Gt), cse=true)
        @test occursin("double u, double J, double c", code)
       


        @syms w T d # Define symbolic variables
        h = 0.01
        G = tf([w^2], [1, 2*d*w, w^2]) * tf(1, [T, 1])
        Gd = tustin(G, h) # Discretize 
        code = SymbolicControlSystems.ccode(Gd, cse=true)
        path = mktempdir()
        filename = joinpath(path, "code.c")
        outname = joinpath(path, "test.so")
        write(joinpath(path, filename), code)
        run(`gcc $filename -lm -shared -o $outname`)

        function c_lsim(u, T, d, w)
            Libc.Libdl.dlopen(outname) do lib
                fn = Libc.Libdl.dlsym(lib, :transfer_function)
                map(u) do u
                    @ccall $(fn)(u::Float64, T::Float64, d::Float64, w::Float64)::Float64
                end
            end
        end

        function c_lsim_ss(u, T, d, w)
            y = zeros(1)
            Libc.Libdl.dlopen(outname) do lib
                fn = Libc.Libdl.dlsym(lib, :transfer_function)
                map(u) do u
                    @ccall $(fn)(y::Ref{Cdouble}, u::Float64, T::Float64, d::Float64, w::Float64)::Cvoid
                    y[]
                end
            end
        end

        function c_lsim_ss_multiinput(u, T, d, w)
            y = zeros(1)
            Y = Libc.Libdl.dlopen(outname) do lib
                fn = Libc.Libdl.dlsym(lib, :transfer_function)
                map(eachcol(u)) do u
                    @ccall $(fn)(y::Ref{Cdouble}, u::Ref{Cdouble}, T::Float64, d::Float64, w::Float64)::Cvoid
                    y[]
                end
            end
            Y'
        end

        function c_lsim_ss(u)
            Y = Libc.Libdl.dlopen(outname) do lib
                fn = Libc.Libdl.dlsym(lib, :transfer_function)
                map(u) do u
                    y = zeros(2)
                    @ccall $(fn)(y::Ref{Cdouble}, u::Float64)::Cvoid
                    y
                end
            end
            reduce(hcat, Y)'
        end

        u = randn(1,1000); # Random input signal 
        T_, d_, w_ = 0.03, 0.2, 2.0 # Define system parameters
        y = c_lsim( u,  T_,  d_,  w_); # Filter u through the C-function filter
        Gd_ = sym2num(Gd, h, Pair.((T, d, w), (T_, d_, w_))...) # Replace symbols with numeric constants
        y_,_ = lsim(ss(Gd_), u); # Filter using Julia
        @test norm(y-y_)/norm(y_) < 1e-10
            

        code = SymbolicControlSystems.ccode(ss(Gd), cse=true)
        write(joinpath(path, filename), code)
        run(`gcc $filename -lm -shared -o $outname`)
        y = c_lsim_ss( u,  T_,  d_,  w_); # Filter u through the 
        y_,_ = lsim(ss(Gd_), u);
        @test norm(y-y_)/norm(y_) < 1e-10 # TODO: figure out why this is more sensitive

        G_ = sym2num(G, Pair.((T, d, w), (T_, d_, w_))...) 
        y_,_ = lsim(c2d(ss(G_), h, :tustin), u);
        @test norm(y-y_)/norm(y_) < 1e-10 # TODO: figure out why this is more sensitive

        ## Multiple inputs
        @syms w T d # Define symbolic variables
        h = 0.01
        G = tf([w^2], [1, 2*d*w, w^2]) * tf(1, [T, 1])
        Gd = tustin(G, h) # Discretize 
        Gd = [Gd 2Gd]
        u = randn(2,1000); # Random input signal 
        T_, d_, w_ = 0.03, 0.2, 2.0 # Define system parameters            

        code = SymbolicControlSystems.ccode(ss(Gd), cse=true)
        write(joinpath(path, filename), code)
        run(`gcc $filename -lm -shared -o $outname`)
        y = c_lsim_ss_multiinput( u,  T_,  d_,  w_); # Filter u through the 
        Gd_ = sym2num(ss(Gd), Pair.((T, d, w), (T_, d_, w_))...) # Replace symbols with numeric constants
        y_,_ = lsim(ss(Gd_), u);
        @test norm(y-y_)/norm(y_) < 1e-10 # TODO: figure out why this is more sensitive


        # test without symbols

        J = 2
        c = 1
        G = tf(1.,[J^2,c,1])
        Gn = tf(1.,[1,1,1])
        Gt = tustin(G, 0.1)
        code = SymbolicControlSystems.ccode(Gt)
        @test occursin("static double u[3] = {0};", code)
        @test occursin("static double y[3] = {0};", code)
        @test occursin("double ui", code)
        
        code = SymbolicControlSystems.ccode(ss(Gt), cse=false)
        @test occursin("double u", code)
       
        code = SymbolicControlSystems.ccode(ss(Gt), cse=true)
        @test occursin("double u", code)

        ## test multi output
        Gt2 = c2d([ss(G); ss(Gn)], 0.1, :tustin)
        code = SymbolicControlSystems.ccode(Gt2, cse=true)
        write(joinpath(path, filename), code)
        run(`gcc $filename -lm -shared -o $outname`)
        @test occursin("double *y, double u", code)

        u = randn(1,1000); # Random input signal 
        y_,_ = lsim(Gt2, u); # Filter using Julia
        y = c_lsim_ss(u);
        @test norm(y'-y_)/norm(y_) < 1e-10

    end




    @testset "Symbolics: Sym -> tf and vice versa" begin
        @info "Symbolics: Testing Sym -> tf and vice versa"

        Symbolics.@variables a b c s z
        @test_both_symb tf([a], [b, c])  a/(b*s + c)
        @test_both_symb tf([a], [b, c], 0.1)  a/(b*z + c)

        @test_both_symb tf([a], [b, c, 1])  a/(b*s^2 + c*s + 1)
        @test_both_symb tf([a], [b, c, 1], 0.1)  a/(b*z^2 + c*z + 1)

        # Test that denominator is still monic
        @test string(to_num(tf([a], [1, c, b]))) == string(a/(s^2 + c*s + b))
        @test string(to_num(tf([a], [1, c, b], 0.1))) == string(a/(z^2 + c*z + b))


        # Test that denominator is still monic
        @test string(to_num(tf([b, a], [1, c, 1]))) == string((b*s + a)/(s^2 + c*s + 1))
        @test string(to_num(tf([b, a], [1, c, 1], 0.1))) == string((b*z + a)/(z^2 + c*z + 1))


        # Test that denominator is still monic
        @test_both_symb tf([a], [1, c, b])  a/(s^2 + c*s + b)
        @test_both_symb tf([a], [1, c, b], 0.1)  a/(z^2 + c*z + b)


        # Test that denominator is still monic
        @test_both_symb tf([b, a], [1, c, 1])  (b*s + a)/(s^2 + c*s + 1)
        @test_both_symb tf([b, a], [1, c, 1], 0.1)  (b*z + a)/(z^2 + c*z + 1)

        # MIMO
        to_num(tf([a], [1, c, b]) .* LinearAlgebra.diagm([1, 1]))
        
    end

    @testset "conversion between SymPy and Symbolics" begin
        @info "Testing conversion between SymPy and Symbolics"

        Gsymbolcs = let
            Symbolics.@variables a b c s
            tf([a], [b, c])
        end
        Gsympy = let
            @syms a b c s
            tf([a], [b, c])
        end

        @test Gsymbolcs == tf(to_num(Gsympy))
        @test Symbolics.symbolics_to_sympy(to_num(Gsymbolcs)) == Sym(Gsympy)
    end
end
