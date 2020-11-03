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
        simplify((C*inv(z*I(size(A,1)) - A)*B + D)[1])
    else
        simplify((C*inv(s*I(size(A,1)) - A)*B + D)[1])
    end
end

function SymPy.Sym(sys::TransferFunction)
    ControlSystems.numeric_type(sys) <: Sym || throw(MethodError(Sym, sys))
    Sym(ss(sys))
end

function ControlSystems.tf(sys::Sym, h=nothing)
    sys = simplify(sys)
    sp.monic(inv(sys)) |> inv |> simplify # makes the denominator monic
    n,d = sp.fraction(sys) .|> sp.Poly
    if h === nothing
        tf(expand_coeffs(n, s), expand_coeffs(d, s))
    else
        tf(expand_coeffs(n, z), expand_coeffs(d, z), h)
    end
end

function expand_coeffs(n, var=s; numeric=false)
    n = sp.Poly(n,var)
    deg = n.degree() |> Float64 |> Int
    c = n.coeffs() # This fails if the coeffs are symbolic
    numeric && (c = Float64.(c))
    [c; zeros(eltype(c), deg-length(c)+1)]
end
expand_coeffs(n::Real, args...; numeric=false) = n


function ControlSystems.minreal(sys::TransferFunction{<:Any, <:ControlSystems.SisoTf{<:Sym}})
    Sym(sys) |> simplify |> tf
end


function sym2num(P, ps, p)
    @warn "TODO: accept a vector of pairs instead"
    # for i = 1:2
    for sym = ps
        ssym = Symbol(sym)
        hasfield(typeof(p), ssym) || continue
        val = float(getfield(p, ssym))
        P = subs(P, (sym, val))
    end
    # end
    n, d = sp.fraction(P) .|> sp.Poly
    nn = expand_coeffs(n, numeric=true)
    nd = expand_coeffs(d, numeric=true)
    G = tf(nn, nd isa Number ? [nd] : nd)
    # isstable(G) || @error "Unstable transfer function with poles $(pole(G))"
    G
end


sym2num(P::TransferFunction, args...) = sym2num(Sym(P), args...)

function sym2num(P::Sym, h::Real, pairs::Pair...)
    # for i = 1:2
    for (sym, val) = pairs
        P = subs(P, (sym, val))
    end
    @show P
    n, d = sp.Poly.(sp.fraction(P), z)
    nn = float.(expand_coeffs(n, z))
    nd = float.(expand_coeffs(d, z))
    G = tf(nn, nd isa Number ? [nd] : nd, h)
    G
end

function sym2num(P::Sym, pairs::Pair...)
    # for i = 1:2
    P.has(z) && error("Found `z` in symbolic expression, provide sample time as second argument to `sym2num`")
    for (sym, val) = pairs
        P = subs(P, (sym, val))
    end
    n, d = sp.Poly.(sp.fraction(P), s)    
    nn = float.(expand_coeffs(n, s))
    nd = float.(expand_coeffs(d, s))
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

function tustin(G::Sym, h)
    G = G.subs(s, (2.0/h) * (z-1)/(z+1))
    Gt = simplify(expand(G))
    tf(Gt, h)
end


function ccode(G::TransferFunction)
    P = Sym(G)
    P.has(z) || error("Did not find `z` in symbolic expression")
    P.has(s) && error("Found `s` in symbolic expression, provide expression in `z`")
    n,d = numvec(G)[], denvec(G)[]
    # @assert d[1] == 1 "Provide monic denominator polynomial"
    if d[1] != 1
        d = simplify.(d ./ d[1])
        n = simplify.(n ./ d[1])
    end
    vars = P.free_symbols
    vars.remove(z)
    vars = collect(vars)
    var_str = ""
    for var in vars
        var_str *= ", double $(var)"
    end

    nu,ny = length(n), length(d)
    code = """
#include <stdio.h>\n
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
    for (i,n) in enumerate(n)
        code *= "    y[0] += $(sp.ccode(n))*u[$(i-1)];\n"
    end
    for (i,n) in enumerate(d[2:end])
        code *= "    y[0] += $(sp.ccode(simplify(-n)))*y[$(i)];\n"
    end
    code *= "    return y[0];\n}"
    print(code)
    clipboard(code)
    code
end


end

