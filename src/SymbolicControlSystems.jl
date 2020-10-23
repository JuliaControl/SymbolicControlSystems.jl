module SymbolicControlSystems
using LinearAlgebra
using ControlSystems, SymPy, Latexify

export sp, sym2num, latextf, Sym, @vars, simplify

const sp = SymPy.PyCall.PyNULL()


function __init__()
    @eval global const s = Sym("s")
    copy!(sp, SymPy.sympy)
end

function SymPy.Sym(sys::StateSpace{<:Any, Sym})
    A,B,C,D = ControlSystems.ssdata(sys)
    simplify((C*inv(s*I(size(A,1)) - A)*B + D)[1])
end

function ControlSystems.tf(sys::Sym)
    sys = simplify(sys)
    sp.monic(inv(sys)) |> inv |> simplify # makes the denominator monic
    n,d = sp.fraction(sys) .|> sp.Poly
    tf(expand_coeffs(n), expand_coeffs(d))
end

function expand_coeffs(n; numeric=false)
    n = sp.Poly(n,s)
    deg = n.degree() |> Float64 |> Int
    c = n.coeffs() # This fails if the coeffs are symbolic
    numeric && (c = Float64.(c))
    [c; zeros(eltype(c), deg-length(c)+1)]
end
expand_coeffs(n::Real; numeric=false) = n


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


function sym2num(P, pairs::Pair...)
    # for i = 1:2
    for (sym, val) = pairs
        ssym = Symbol(sym)
        P = subs(P, (sym, val))
    end
    n, d = sp.fraction(P) .|> sp.Poly
    nn = expand_coeffs(n)
    nd = expand_coeffs(d)
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


end
