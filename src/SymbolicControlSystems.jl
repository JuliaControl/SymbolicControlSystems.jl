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
    if isdiscrete(sys)
        ssys = sp.Poly(numvec(sys)[], z)/sp.Poly(denvec(sys)[], z)
    else
        ssys = sp.Poly(numvec(sys)[], s)/sp.Poly(denvec(sys)[], s)
    end
end

"""
    ControlSystems.tf(sys::Sym, h = nothing)

Convert a symbolic, rational expression into a transfer function. `h` denotes the sample time if the system is discrete.
"""
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

function ControlSystems.tf(sys::StateSpace{<:Any, Sym})
    n,p = simplify.(sp.Poly.(simplify.(sp.fraction(simplify(Sym(sys)))), s))
    tf(simplify(n/p))
end


sym2num(P::TransferFunction, args...) = sym2num(Sym(P), args...)

"""
    sym2num(G, h::Real, pairs::Pair...)

Replace symbols by numbers. `h` indicates the sample rate if the system is discrete.
"""
function sym2num(P::Sym, h::Real, pairs::Pair...)
    P.has(s) && error("Found `z` in symbolic expression, provide sample time as second argument to `sym2num`")
    for (sym, val) = pairs
        P = subs(P, (sym, val))
    end
    n, d = sp.Poly.(sp.fraction(P), z)
    nn = float.(expand_coeffs(n, z))
    nd = float.(expand_coeffs(d, z))
    if nd[1] != 1 # make monic
        nn,nd = monic(nn,nd)
    end
    G = tf(nn, nd isa Number ? [nd] : nd, h)
    G
end


"""
    sym2num(G, pairs::Pair...)

Replace symbols by numbers.
"""
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
        nn,nd = monic(nn,nd)
    end
    G = tf(nn, nd isa Number ? [nd] : nd)
    G
end


"""
    tustin(G::LTISystem, h)

Discretize `G` using the Tustin (bi-linear) method for sample time `h`.
"""
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
    n,d = monic(n,d)
    tf(n,d,h)
end

"""
    doubleeuler(sys::AbstractStateSpace{<:Continuous}, Ts)

Discretize `sys` with a second-order approximation to the exponential map (Forward Euler is a first-order approximation). This is useful if the symbolic expressions for the true ZoH-discretization becomes too complicated. A second-order approximation is in many cases indistinguishable from the true ZoH-discretization.

# Example
```julia
w = 2pi .* exp10.(LinRange(-1, log10(25), 400))
sys = ssrand(1,1,4)
bodeplot(sys, w, lab="cont")
bodeplot!(c2d(sys, 0.001, :fwdeuler), w, label="fwdeuler")
bodeplot!(doubleeuler(sys, 0.001), w, label="double euler")
```
"""
function doubleeuler(sys::AbstractStateSpace{<:Continuous}, Ts)
    A, B, C, D = ssdata(sys)
    T = promote_type(eltype.((A,B,C,D))...)
    ny, nu = size(sys)
    nx = sys.nx
    Ad, Bd, Cd, Dd = (I + Ts*A + Ts^2/2*A^2), Ts*B, C, D
    ss(Ad, Bd, Cd, Dd, Ts)
end


"""
    ccode(G; simplify = identity, cse = true)

Return a string with C-code for filtering a signal `u` through `G`. 

# Arguments:
- `G`: A linear system
- `simplify`: A function for symbolic simplification. You may try `Sympy.simplify`, but for large systems, this will take a long time to compute.
- `cse`: Perform common subexpression elimination. This generally improvems the performance of the generated code.
"""
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
        if length(final) == 1
            final = final[]
        end
        n = final[1:length(n)]
        d = final[length(n)+1:end]
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
    println(code)
    clipboard(code)
    code
end

function ccode(sys::StateSpace{<:Discrete}; cse=true, function_name = "transfer_function")
    sys.nu == 1 || throw(ArgumentError("Multiple input not yet supported"))
    nx = sys.nx
    ny = sys.ny
    u = Sym("u")
    x = [Sym("x[$(i-1)]") for i in 1:nx]
    # @show P
    if ControlSystems.numeric_type(sys) <: SymPy.Sym
        P = Sym(sys)
        vars = P.free_symbols
        vars.remove(z)
        vars = collect(vars)
        vars = sort(vars, by=string)
        var_str = ""
        for var in vars
            var_str *= ", double $(var)"
        end
    else
        var_str = ""
    end
    x1 = zeros(Sym, nx) # workaround for strange bug with undefinied referece appearing in Pluto only
    y = zeros(Sym, ny)
    @show x1 = mul!(x1, sys.A, x) + sys.B*u
    @show y = mul!(y, sys.C, x) + sys.D*u
    # @show y = sp.collect.(y, x)
    
    code = """
#include <stdio.h>\n
#include <math.h>\n
void $(function_name)(double *y, double u$(var_str)) {
    static double x[$(nx)] = {0};  // Current state
    double xp[$(nx)] = {0};        // Next state
    int i;
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
        code *= "    y[$(i-1)] = ($(sp.ccode(n)));\n"
    end
    code *= """

        // Make the predicted state the current state
        for (i=0; i < $(nx); ++i) {
            x[i] = xp[i];
        }
    """
    code *= "\n}"

    println(code)
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
    println(code)
    clipboard(code)
    code
end


latextf(x::TransferFunction, args...) = latextf(Sym(x), args...)

Latexify.Latexify(x::TransferFunction) = latextf(x, args...)

function monic(n,d) 
    n = n ./ d[1]
    d = d ./ d[1] 
    if n isa Number
        n = [n]
    end
    n,d
end

"""
    latextf(x, mon = true)

Return a latex string representing `x`. If `mon`, then the denominator polynomial will be made monic (leading coefficient = 1).
"""
function latextf(x::Sym, mon=true)
    var = x.has(z) ? z : s
    n,d = sp.fraction(x)
    if mon
        n,d = monic(n,d)
    end
    n,d = sp.collect.((n,d), var)
    str = "\$\\dfrac{$((n))}{$((d))}\$"
    str = replace(str, '*' => "")
    clipboard(str)
    str
end

show_construction(sys::LTISystem) = show_construction(stdout, sys)
function show_construction(io::IO, sys::LTISystem)
    sys = StateSpace(sys)
    println(io, "A = ", sys.A)
    println(io, "B = ", sys.B)
    println(io, "C = ", sys.C)
    println(io, "D = ", sys.D)
    if isdiscrete(sys)
        println(io, "sys = ss(A,B,C,D,$(sys.Ts))")
    else
        println(io, "sys = ss(A,B,C,D)")
    end
end

"""
sys2vec = @(sys) [
        size(sys.A,1)
        size(sys.B,2)
        size(sys.C,1)
        sys.A(:)
        sys.B(:)
        sys.C(:)
        sys.D(:)
    ]
"""
function vec2sys(v, ts=nothing)
    nx = Int(v[1])
    nu = Int(v[2])
    ny = Int(v[3])
    ai = (1:nx^2) .+ 3
    bi = (1:nx*nu) .+ ai[end]
    ci = (1:nx*ny) .+ bi[end]
    di = (1:nu*ny) .+ ci[end]
    A = reshape(v[ai], nx, nx)
    B = reshape(v[bi], nx, nu)
    C = reshape(v[ci], ny, nx)
    D = reshape(v[di], ny, nu)
    ts === nothing ? ss(A,B,C,D) : ss(A,B,C,D, ts)
end


print_c_array(a::AbstractArray, args...; kwargs...) = print_c_array(stdout, a, args...; kwargs...)
function print_c_array(io, a::AbstractVector, name="vec"; cse=false)
    l = length(a)
    if cse
        a = write_array_cse(io, a, name)
    end
    println(io, "    double $name[$l];")
    for i = 1:l
        println(io, "    $name[$(l-1)] = $(a[i]);")
    end
end
function print_c_array(io, a::AbstractMatrix, name="mat"; cse=false)
    r,c = size(a)
    if cse
        a = write_array_cse(io, a, name)
    end
    println(io, "    double $name[$r][$c];")
    for i = 1:r, j = 1:c
        println(io, "    $name[$(i-1)][$(j-1)] = $(a[i,j]);")
    end
end

function print_c_array_interp(io, a::AbstractMatrix, t, name="mat"; cse=false)
    r,c = size(a)
    if cse
        a = write_array_cse(io, a, name)
    end
    println(io, "    double $name[$r][$c];")
    for i = 1:r, j = 1:c
        println(io, "    $name[$(i-1)][$(j-1)] = $(a[i,j]);")
    end
end

function write_array_cse(io, a, name="x")
    subex, final = sp.cse(a, symbols=[SymPy.Sym(name*string(i)) for i in 1:200])
    new_a = final[]
    for se in subex
        println(io, "    double $(se[1]) = $(sp.ccode(se[2]));")
    end
    new_a
end

end

