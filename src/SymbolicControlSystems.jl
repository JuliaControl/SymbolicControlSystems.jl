__precompile__(false)
module SymbolicControlSystems
using LinearAlgebra
using ControlSystems, SymPy, Latexify
import Symbolics as Symb
import Symbolics: Num
using InteractiveUtils

export sp,
    sym2num,
    latextf,
    Sym,
    @vars,
    simplify,
    tustin,
    doubleeuler,
    s,
    z,
    print_c_array,
    show_construction

# export symbolics2sym

const sp = SymPy.PyCall.PyNULL()
const s = SymPy.Sym("s")
const z = SymPy.Sym("z")

const NumOrDiv = Union{Num, Symb.SymbolicUtils.Div}


function __init__()
    # global const s = Sym("s")
    # global const z = Sym("z")
    copy!(sp, SymPy.sympy)
end


function SymPy.Sym(sys::StateSpace{<:Any,Sym})
    A, B, C, D = ControlSystems.ssdata(sys)
    if isdiscrete(sys)
        (C*inv(z * I(size(A, 1)) - A)*B+D)[1]
    else
        (C*inv(s * I(size(A, 1)) - A)*B+D)[1]
    end
end

function SymPy.Sym(sys::TransferFunction)
    if isdiscrete(sys)
        ssys = sp.Poly(numvec(sys)[], z) / sp.Poly(denvec(sys)[], z)
    else
        ssys = sp.Poly(numvec(sys)[], s) / sp.Poly(denvec(sys)[], s)
    end
end

function Num(sys::StateSpace{<:Any,Num})
    A, B, C, D = ControlSystems.ssdata(sys)
    λ = isdiscrete(sys) ? Symb.@variables(z) : Symb.@variables(s)
    λ = λ[]
    Symb.simplify((C*inv(λ * I(size(A, 1)) - A)*B+D)[1])
end

function Num(sys::TransferFunction)
    λ = isdiscrete(sys) ? Symb.@variables(z) : Symb.@variables(s)
    λ = λ[]
    num = sum(((i, t),) -> t * λ^(i-1), enumerate(reverse(numvec(sys)[]))) |> Symb.simplify
    den = sum(((i, t),) -> t * λ^(i-1), enumerate(reverse(denvec(sys)[]))) |> Symb.simplify
    Symb.simplify(num/den)
end

"""
    ControlSystems.tf(sys::Sym, h = nothing)

Convert a symbolic, rational expression into a transfer function. `h` denotes the sample time if the system is discrete.
"""
function ControlSystems.tf(sys::Sym, h = nothing)
    n, d = sp.fraction(sys)
    d = sp.Poly(d)
    # d = d.monic() # Don't do this here
    n = n isa Number ? n : sp.Poly(n)
    if h === nothing
        tf(expand_coeffs(n, s), expand_coeffs(d, s))
    else
        tf(expand_coeffs(n, z), expand_coeffs(d, z), h)
    end
end

function expand_coeffs(n, var; numeric = false)
    n = sp.Poly(n, var)
    deg = n.degree() |> Float64 |> Int
    c = n.all_coeffs() # This fails if the coeffs are symbolic
    numeric && (c = Float64.(c))
    [c; zeros(eltype(c), deg - length(c) + 1)]
end
expand_coeffs(n::Real, args...; numeric = false) = n

function ControlSystems.tf(sys::NumOrDiv, h = nothing)
    sp = Symb.symbolics_to_sympy(sys)
    G = tf(sp, h)
    tf(Num.(numvec(G)[]), Num.(denvec(G)[]), G.timeevol)
end

Base.:(==)(s1::TransferFunction{<:Any,<:ControlSystems.SisoTf{Num}}, s2::TransferFunction{<:Any,<:ControlSystems.SisoTf{<:NumOrDiv}}) = isequal(Num(s1), Num(s2))

Base.promote_op(::typeof(/),::Type{NumOrDiv},::Type{NumOrDiv}) = Num # This is required to make conversion to ss work. Arithmetic operaitons on Num are super type unstable so inference fails https://github.com/JuliaSymbolics/Symbolics.jl/issues/626

function ControlSystems.minreal(sys::TransferFunction{<:Any,<:ControlSystems.SisoTf{<:Sym}})
    Sym(sys) |> simplify |> tf
end

function ControlSystems.tf(sys::StateSpace{<:Any,Sym})
    n, p = simplify.(sp.Poly.(simplify.(sp.fraction(simplify(Sym(sys)))), s))
    tf(simplify(n / p))
end


function ControlSystems.minreal(sys::StateSpace{<:Any,NumOrDiv})
    # sys |> Symb.Num .|> Symb.symbolics_to_sympy .|> sp.simplify
    nsys = Num(sys)
    nsys = Symb.simplify.(nsys)
    nsys = Symb.simplify_fractions.(nsys)
end

function Num(x::Sym)
    try
        return Float64(x)
    catch
        Symb.Num(Symb.variable(Symbol(x); T=Real))
    end
end


sym2num(P::TransferFunction, args...) = sym2num(Sym(P), args...)

"""
    sym2num(G, h::Real, pairs::Pair...)

Replace symbols by numbers. `h` indicates the sample rate if the system is discrete.
"""
function sym2num(P::Sym, h::Real, pairs::Pair...)
    P.has(s) && error(
        "Found `z` in symbolic expression, provide sample time as second argument to `sym2num`",
    )
    for (sym, val) in pairs
        P = subs(P, (sym, val))
    end
    n, d = sp.Poly.(sp.fraction(P), z)
    nn = float.(expand_coeffs(n, z))
    nd = float.(expand_coeffs(d, z))
    if nd[1] != 1 # make monic
        nn, nd = monic(nn, nd)
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
    P.has(z) && error(
        "Found `z` in symbolic expression, provide sample time as second argument to `sym2num`",
    )
    P.has(s) || error("Found no `s` in symbolic expression, provide")
    for (sym, val) in pairs
        P = subs(P, (sym, val))
    end
    n, d = sp.Poly.(sp.fraction(P), s)
    nn = float.(expand_coeffs(n, s))
    nd = float.(expand_coeffs(d, s))
    if nd[1] != 1 # make monic
        nn, nd = monic(nn, nd)
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
    Gt = Gs.subs(s, hh * (z - 1) / (z + 1)) |> simplify
    @info "Perforing pole-zero cancellation (this might take a while)"
    @time Gt = sp.cancel(Gt) # Simplify takes forever with floating point values
    @time Gt = sp.collect(Gt, z) # Simplify takes forever with floating point values
    Gt = Gt.subs(hh, (2 / h))
    Gtf = tf(Gt, h)
    n, d = numvec(Gtf)[], denvec(Gtf)[]
    d[1] == 1 && (return Gtf)
    n, d = monic(n, d)
    tf(n, d, h)
end

"""
    doubleeuler(sys::AbstractStateSpace{<:Continuous}, Ts)

Discretize `sys` with a second-order approximation to the exponential map (Forward Euler is a first-order approximation). This is useful if the symbolic expressions for the true ZoH-discretization becomes too complicated. A second-order approximation is in many cases indistinguishable from the true ZoH-discretization, but requires a well-balanced state-space realization.

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
    T = promote_type(eltype.((A, B, C, D))...)
    ny, nu = size(sys)
    nx = sys.nx
    A = [A B; zeros(nu, nx + nu)]
    A2 = A*A
    # A3 = A2*A
    M = I + Ts * A + Ts^2 / 2 * A2 #+ Ts^3 / 6 * A3
    Ad = M[1:nx, 1:nx]
    Bd = M[1:nx, nx+1:nx+nu]
    ss(Ad, Bd, C, D, Ts)
end


"""
    ccode(G; simplify = identity, cse = true)

Return a string with C-code for filtering a signal `u` through `G`. 

# Arguments:
- `G`: A linear system
- `simplify`: A function for symbolic simplification. You may try `Sympy.simplify`, but for large systems, this will take a long time to compute.
- `cse`: Perform common subexpression elimination. This generally improvems the performance of the generated code.
"""
function ccode(G::TransferFunction; simplify = identity, cse = true)
    P = Sym(G)
    P.has(z) || error("Did not find `z` in symbolic expression")
    P.has(s) && error("Found `s` in symbolic expression, provide expression in `z`")
    n, d = numvec(G)[], denvec(G)[]
    if d[1] != 1
        @info "Denominator polynomial not monic, dividing by the leading coefficient."
        n = simplify.(n ./ d[1])
        d = simplify.(d ./ d[1])
    end
    @info "Calling free_symbols"
    @info vars = P.free_symbols
    vars.remove(z)
    vars = collect(vars)
    vars = sort(vars, by = string)
    var_str = ""
    for var in vars
        var_str *= ", double $(var)"
    end

    nu, ny = length(n), length(d)
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
        subex, final = sp.cse([n; d])
        if length(final) == 1
            final = final[]
        end
        n = final[1:length(n)]
        d = final[length(n)+1:end]
        for se in subex
            code *= "    double $(se[1]) = $(sp.ccode(se[2]));\n"
        end
    end
    for (i, n) in enumerate(n)
        @info "Writing numerator term $i/$(length(n))"
        code *= "    y[0] += ($(sp.ccode(n)))*u[$(i-1)];\n"
    end
    for (i, n) in enumerate(d[2:end])
        @info "Writing denominator term $i/$(length(d)-1)"
        code *= "    y[0] += ($(sp.ccode(-n)))*y[$(i)];\n"
    end
    code *= "    return y[0];\n}"
    println(code)
    try
        clipboard(code)
    catch
    end
    code
end

function ccode(sys::StateSpace{<:Discrete}; cse = true, function_name = "transfer_function")
    sys.nu == 1 || throw(ArgumentError("Multiple input not yet supported"))
    nx = sys.nx
    ny = sys.ny
    u = Sym("u")
    x = [Sym("x[$(i-1)]") for i = 1:nx]
    # @show P
    if ControlSystems.numeric_type(sys) <: SymPy.Sym
        P = Sym(sys)
        vars = P.free_symbols
        vars.remove(z)
        vars = collect(vars)
        vars = sort(vars, by = string)
        var_str = ""
        for var in vars
            var_str *= ", double $(var)"
        end
    else
        var_str = ""
    end
    x1 = zeros(Sym, nx) # workaround for strange bug with undefinied referece appearing in Pluto only
    y = zeros(Sym, ny)
    @show x1 = mul!(x1, sys.A, x) + sys.B * u
    @show y = mul!(y, sys.C, x) + sys.D * u
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
        subex, final = sp.cse([x1; y])
        x1 = final[][1:length(x1)]
        y = final[][length(x1)+1:end]
        for se in subex
            code *= "    double $(se[1]) = $(sp.ccode(se[2]));\n"
        end
    end
    code *= "\n    // Advance the state xp = Ax + Bu\n"
    for (i, n) in enumerate(x1)
        code *= "    xp[$(i-1)] = ($(sp.ccode(n)));\n"
    end
    code *= """

        // Accumulate the output y = C*x + D*u
    """

    for (i, n) in enumerate(y)
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
    try
        clipboard(code)
    catch
    end
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
    try
        clipboard(code)
    catch
    end
    code
end


latextf(x::TransferFunction, args...) = latextf(Sym(x), args...)

Latexify.Latexify(x::TransferFunction) = latextf(x, args...)

function monic(n, d)
    n = n ./ d[1]
    d = d ./ d[1]
    if n isa Number
        n = [n]
    end
    n, d
end

"""
    latextf(x, mon = true)

Return a latex string representing `x`. If `mon`, then the denominator polynomial will be made monic (leading coefficient = 1).
"""
function latextf(x::Sym, mon = true)
    var = x.has(z) ? z : s
    n, d = sp.fraction(x)
    if mon
        n, d = monic(n, d)
    end
    n, d = sp.collect.((n, d), var)
    str = "\$\\dfrac{$((n))}{$((d))}\$"
    str = replace(str, '*' => "")
    try
        clipboard(code)
    catch
    end
    str
end

show_construction(sys::LTISystem; kwargs...) = show_construction(stdout, sys; kwargs...)
function show_construction(io::IO, sys::LTISystem; letb = true)
    # sys = StateSpace(sys)
    letb && println(io, "sys = let")
    prestr = letb ? "    " : "" 
    println(io, prestr*"tempA = ", sys.A)
    println(io, prestr*"tempB = ", sys.B)
    println(io, prestr*"tempC = ", sys.C)
    println(io, prestr*"tempD = ", sys.D)
    letb || print(io, "sys = ")
    if isdiscrete(sys)
        println(io, prestr*"ss(tempA, tempB, tempC, tempD, $(sys.Ts))")
    else
        println(io, prestr*"ss(tempA, tempB, tempC, tempD)")
    end
    letb && println(io, "end")
    nothing
end

"""
sys2vec = @(sys) [
        size(sys.A,1);
        size(sys.B,2);
        size(sys.C,1);
        sys.A(:);
        sys.B(:);
        sys.C(:);
        sys.D(:)
    ]
"""
function vec2sys(v, ts = nothing)
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
    ts === nothing ? ss(A, B, C, D) : ss(A, B, C, D, ts)
end


print_c_array(a::AbstractArray, args...; kwargs...) =
    print_c_array(stdout, a, args...; kwargs...)
function print_c_array(io, a::AbstractVector, name = "vec"; cse = false, s = "", struct_name=nothing)
    l = length(a)
    a = float.(a)
    if cse
        a = write_array_cse(io, a, name, s)
    end
    if struct_name === nothing
        s == "" && println(io, "double $name[$l];")
        struct_name = ""
    else
        struct_name = struct_name*"->"
    end
    for i = 1:l
        println(io, "$(struct_name)" * s * "$name[$(i-1)] = $(a[i]);")
    end
    println(io)
end
function print_c_array(io, a::AbstractMatrix, name = "mat"; cse = false, s = "", struct_name=nothing)
    r, c = size(a)
    a = float.(a)
    if cse
        a = write_array_cse(io, a, name, s)
    end
    if struct_name === nothing
        s == "" && println(io, "double $name[$r][$c];")
        struct_name = ""
    else
        struct_name = struct_name*"->"
    end
    for i = 1:r, j = 1:c
        println(io, "$(struct_name)" * s * "$name[$(i-1)][$(j-1)] = $(a[i,j]);")
    end
    println(io)
end

function print_c_array(io, a::AbstractArray{<:Any,3}, name = "array"; cse = false, s = "", struct_name=nothing)
    r, c, d = size(a)
    a = float.(a)
    if cse
        a = write_array_cse(io, a, name, s)
    end
    if struct_name === nothing
        s == "" && println(io, "double $name[$r][$c][$d];")
        struct_name = ""
    else
        struct_name = struct_name*"->"
    end
    for i = 1:r, j = 1:c, k = 1:d
        println(io, "$(struct_name)" * s * "$name[$(i-1)][$(j-1)][$(k-1)] = $(a[i,j,k]);")
    end
    println(io)
end


# Interpolation
function print_c_array(
    io,
    a::Vector{<:AbstractArray},
    t::AbstractVector,
    name          = "mat";
    cse           = false,
    s             = "",
    print_vector  = true,
    print_logic   = true,
    struct_name::Union{Nothing, String} = nothing,
    struct_type   = nothing,
    ivecname = name * "_interp_vect",
)
    length(a) == length(t) || throw(ArgumentError("length mismatch between a and t"))
    print_vector && print_c_array(io, t, ivecname; cse, s, struct_name)
    # for i in eachindex(t)
    #     iname = name*"_$(i-1)"
    #     print_c_array(io, a[i], iname; cse, s)
    # end
    a = cat(a..., dims = 3)

    if struct_name !== nothing
        struct_type isa String || throw(ArgumentError("Expected a struct type since you provided struct name"))
        print_c_array(io, a, name * "_interp"; cse, s, struct_name)
        println(io, "void interpolate_$(name)($(struct_type) *$(struct_name), double t) {")
        ivecname = struct_name*"->"*ivecname
    else
        print_c_array(io, a, name * "_interp"; cse, s)
        println(
        io,
        "void interpolate_$(name)(double *$(name), double *$(name)_interp, double *$(ivecname), double t) {"
        )
    end

    println(
        io,
        """
    int k = 0;
    for (k = 0; k < $(length(t)-2); ++k) { // Loop to find correct interval
        if (t < $(ivecname)[k+1]) {
            break;
        }
    } // t is now between indices k and k+1 modolu edge cases below
    double t0 = $(ivecname)[k];
    double t1 = $(ivecname)[k+1];
    double l = t1 - t0;      // Length of interval between points
    double alpha = (t-t0)/l; // Between 0 and 1, proportion of k+1 point
    if (t < $(ivecname)[0]) { // edge cases
        alpha = 0;
    }else if (t > $(ivecname)[$(length(t)-1)]) {
        alpha = 1;
    }
$(interpolator_string(a, name, struct_name, struct_type))}
""",
    )
end

function interpolator_string(a::AbstractArray{<:Any,3}, name, struct_name, struct_type)
    if struct_name === nothing
        """
            for (int i = 0; i < $(size(a, 1)); ++i) {     // $name has dimensions $((size(a)[1:2]))
                for (int j = 0; j < $(size(a, 2)); ++j) { // $(name)_interp has dimensions $((size(a)))
                    $(name)[i][j] = (1-alpha)*$(name)_interp[i][j][k] + alpha*$(name)_interp[i][j][k+1];
                }
            }
        """
    else
        name = struct_name*"->"*name
        """
            for (int i = 0; i < $(size(a, 1)); ++i) {     // $name has dimensions $((size(a)[1:2]))
                for (int j = 0; j < $(size(a, 2)); ++j) { // $(name)_interp has dimensions $((size(a)))
                    $(name)[i][j] = (1-alpha)*$(name)_interp[i][j][k] + alpha*$(name)_interp[i][j][k+1];
                }
            }
        """
    end
end

# StateSpace
function print_c_array(io, sys::AbstractStateSpace, args...; s = "", en = "", struct_name=nothing)
    print_c_array(io, sys.A, "A" * en; s, struct_name)
    print_c_array(io, sys.B, "B" * en; s, struct_name)
    print_c_array(io, sys.C, "C" * en; s, struct_name)
    print_c_array(io, sys.D, "D" * en; s, struct_name)
end


# Interpolated StateSpace
function print_c_array(
    io,
    sys::Vector{<:AbstractStateSpace},
    t::AbstractVector,
    name = "sys";
    cse = false,
    s = "",
    en = "",
    struct_name::Union{Nothing, String} = nothing,
    struct_type   = nothing,
)
    length(sys) == length(t) || throw(ArgumentError("length mismatch between a and t"))
    ivecname = name * en * "_interp_vect"
    print_c_array(io, t, ivecname; cse, s, struct_name)
    for p in (:A, :B, :C, :D)
        A = getproperty.(sys, p)
        print_c_array(io, A, t, name * "_" * string(p) * en; cse, s, print_vector = false, struct_name, struct_type, ivecname)
    end
end

function write_array_cse(io, a, name = "x", s = "")
    subex, final = sp.cse(a, symbols = [SymPy.Sym(name * string(i)) for i = 1:200])
    new_a = final[]
    for se in subex
        println(io, "double $(se[1]) = $(sp.ccode(se[2]));")
    end
    new_a
end



end