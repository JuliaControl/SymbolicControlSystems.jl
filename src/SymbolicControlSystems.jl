module SymbolicControlSystems
using LinearAlgebra
using ControlSystems, SymPy, Latexify
using InteractiveUtils

export sp, sym2num, latextf, Sym, @vars, simplify, tustin, s, z

const sp = SymPy.PyCall.PyNULL()


function __init__()
    @eval global const s = Sym("s")
    @eval global const z = Sym("z")
    copy!(sp, SymPy.sympy)
end


function SymPy.Sym(sys::StateSpace{<:Any, Sym})
    A,B,C,D = ControlSystems.ssdata(sys)
    if isdiscrete(sys)
        (C*inv(z*I(size(A,1)) - A)*B + D)[1]
    else
        (C*inv(s*I(size(A,1)) - A)*B + D)[1]
    end
end

function SymPy.Sym(sys::TransferFunction)
    ControlSystems.numeric_type(sys) <: Sym || throw(MethodError(Sym, sys))
    if isdiscrete(sys)
        ssys = sp.Poly(numvec(sys)[], z)/sp.Poly(denvec(sys)[], z)
    else
        ssys = sp.Poly(numvec(sys)[], s)/sp.Poly(denvec(sys)[], s)
    end
end

function ControlSystems.tf(sys::Sym, h=nothing)
    n,d = sp.fraction(sys)
    d = sp.Poly(d)
    # d = d.monic() # Don't do this here
    n = n isa Number ? n : sp.Poly(n)
    if h === nothing
        tf(expand_coeffs(n, s), expand_coeffs(d, s))
    else
        tf(expand_coeffs(n, z), expand_coeffs(d, z), h)
    end
end

function expand_coeffs(n, var; numeric=false)
    n = sp.Poly(n,var)
    deg = n.degree() |> Float64 |> Int
    c = n.all_coeffs() # This fails if the coeffs are symbolic
    numeric && (c = Float64.(c))
    [c; zeros(eltype(c), deg-length(c)+1)]
end
expand_coeffs(n::Real, args...; numeric=false) = n


function ControlSystems.minreal(sys::TransferFunction{<:Any, <:ControlSystems.SisoTf{<:Sym}})
    Sym(sys) |> simplify |> tf
end


# function sym2num(P, ps, p)
#     @warn "TODO: accept a vector of pairs instead"
#     # for i = 1:2
#     for sym = ps
#         ssym = Symbol(sym)
#         hasfield(typeof(p), ssym) || continue
#         val = float(getfield(p, ssym))
#         P = subs(P, (sym, val))
#     end
#     # end
#     n, d = sp.fraction(P) .|> sp.Poly
#     nn = expand_coeffs(n, numeric=true)
#     nd = expand_coeffs(d, numeric=true)
#     G = tf(nn, nd isa Number ? [nd] : nd)
#     # isstable(G) || @error "Unstable transfer function with poles $(pole(G))"
#     G
# end


sym2num(P::TransferFunction, args...) = sym2num(Sym(P), args...)

function sym2num(P::Sym, h::Real, pairs::Pair...)
    P.has(s) && error("Found `z` in symbolic expression, provide sample time as second argument to `sym2num`")
    for (sym, val) = pairs
        P = subs(P, (sym, val))
    end
    n, d = sp.Poly.(sp.fraction(P), z)
    nn = float.(expand_coeffs(n, z))
    nd = float.(expand_coeffs(d, z))
    if nd[1] != 1 # make monic
        nn ./= nd[1]
        nd ./= nd[1]
    end
    G = tf(nn, nd isa Number ? [nd] : nd, h)
    G
end

function sym2num(P::Sym, pairs::Pair...)
    # for i = 1:2
    P.has(z) && error("Found `z` in symbolic expression, provide sample time as second argument to `sym2num`")
    P.has(s) || error("Found no `s` in symbolic expression, provide")
    for (sym, val) = pairs
        P = subs(P, (sym, val))
    end
    n, d = sp.Poly.(sp.fraction(P), s)    
    nn = float.(expand_coeffs(n, s))
    nd = float.(expand_coeffs(d, s))
    if nd[1] != 1 # make monic
        nn ./ nd[1]
        nd ./ nd[1]
    end
    G = tf(nn, nd isa Number ? [nd] : nd)
    G
end

function latextf(x)
    @warn "TODO: group coefficients of powers of s together"
    n,d = sp.fraction(x)
    s = "\$\\dfrac{$((n))}{$((d))}\$"
    s = replace(s, '*' => "")
    clipboard(s)
    s
end

# replace_controllers(ex) = subs(ex, (C, Cᵥ2*(s+Cₚ2)), (Cₚ, Cₚ2), (Cᵥ,Cᵥ2)) |> simplify
# numeric_tf(ex) = sym2num(replace_controllers(ex), ps, pn)


function tustin(G::LTISystem, h)
    iscontinuous(G) || throw(DomainError("Cannot discretize a discrete-domain system."))
    tustin(Sym(G), h)
end

function tustin(Gs::Sym, h)
    hh = Sym("hh")
    Gt = Gs.subs(s, hh * (z-1)/(z+1)) |> simplify
    @info "Perforing pole-zero cancellation (this might take a while)"
    @time Gt = sp.cancel(Gt) # Simplify takes forever with floating point values
    @time Gt = sp.collect(Gt, z) # Simplify takes forever with floating point values
    Gt = Gt.subs(hh, (2/h))
    Gtf = tf(Gt, h)
    n,d = numvec(Gtf)[], denvec(Gtf)[]
    d[1] == 1 && (return Gtf)
    n = n ./ d[1]
    d = d ./ d[1]
    tf(n,d,h)
end


function ccode(G::TransferFunction; simplify = identity, cse=true)
    P = Sym(G)
    P.has(z) || error("Did not find `z` in symbolic expression")
    P.has(s) && error("Found `s` in symbolic expression, provide expression in `z`")
    n,d = numvec(G)[], denvec(G)[]
    if d[1] != 1
        @info "Denominator polynomial not monic, dividing by the leading coefficient."
        n = simplify.(n ./ d[1])
        d = simplify.(d ./ d[1])
    end
    @info "Calling free_symbols"
    @info vars = P.free_symbols
    vars.remove(z)
    vars = collect(vars)
    vars = sort(vars, by=string)
    var_str = ""
    for var in vars
        var_str *= ", double $(var)"
    end

    nu,ny = length(n), length(d)
    code = """
#include <stdio.h>\n
#include <math.h>\n
double transfer_function(double ui$(var_str)) {
    static double u[$(nu)] = {0};
    static double y[$(ny)] = {0};
    int i;
    for (i=$(nu-1); i > 0; --i) {
        u[i] = u[i-1];
    }
    u[0] = ui;
    for (i=$(ny-1); i > 0; --i) {
        y[i] = y[i-1];
    }
    y[0] = 0;
"""
    if cse
    @info "Finding common subexpressions"
        subex, final = sp.cse([n;d])
        n = final[][1:length(n)]
        d = final[][length(n)+1:end]
        for se in subex
            code *= "    double $(se[1]) = $(sp.ccode(se[2]));\n"
        end
    end
    for (i,n) in enumerate(n)
        @info "Writing numerator term $i/$(length(n))"
        code *= "    y[0] += ($(sp.ccode(n)))*u[$(i-1)];\n"
    end
    for (i,n) in enumerate(d[2:end])
        @info "Writing denominator term $i/$(length(d)-1)"
        code *= "    y[0] += ($(sp.ccode(-n)))*y[$(i)];\n"
    end
    code *= "    return y[0];\n}"
    print(code)
    clipboard(code)
    code
end

function ccode(sys::StateSpace{<:Discrete}; cse=true)
    nx = sys.nx
    u = Sym("u")
    x = [Sym("x[$(i-1)]") for i in 1:nx]
    P = Sym(sys)
    vars = P.free_symbols
    vars.remove(z)
    vars = collect(vars)
    vars = sort(vars, by=string)
    var_str = ""
    for var in vars
        var_str *= ", double $(var)"
    end
    x1 = sp.collect.(sys.A*x + sys.B*u, z)
    @show y = (sys.C*x + sys.D*u)[]
    @show y = sp.collect(y, z)
    
    code = """
#include <stdio.h>\n
#include <math.h>\n
double transfer_function(double u$(var_str)) {
    static double x[$(nx)] = {0};  // Current state
    double xp[$(nx)] = {0};        // Next state
    int i;
    double y = 0; // Output
"""
    if cse
        @info "Finding common subexpressions"
        subex, final = sp.cse([x1;y])
        x1 = final[][1:length(x1)]
        y = final[][length(x1)+1:end]
        for se in subex
            code *= "    double $(se[1]) = $(sp.ccode(se[2]));\n"
        end
    end
    code *= "\n    // Advance the state xp = Ax + Bu\n"
    for (i,n) in enumerate(x1)
        code *= "    xp[$(i-1)] = ($(sp.ccode(n)));\n"
    end
    code *= """

        // Accumulate the output y = C*x + D*u
    """

    for (i,n) in enumerate(y)
        code *= "    y += ($(sp.ccode(n)));\n"
    end
    code *= """

        // Make the predicted state the current state
        for (i=0; i < $(nx); ++i) {
            x[i] = xp[i];
        }
    """
    code *= "    return y;\n}"

    print(code)
    clipboard(code)
    code
end


function structured_text(code)

    code = replace(code, "pow" => "EXPT")
    code = replace(code, "static " => "")
    
    code = replace(code, r"double (\w+) =" => s"\1 : LREAL :=")
    code = replace(code, r"double (\w+)" => s"\1 : LREAL")
    code = replace(code, r"(.{1,20}?) \+= " => s"\1 := \1 + ")
    code = replace(code, " = " => " := ")
    clipboard(code)
    code

end

end

