# SymbolicControlSystems


Utilities for
- Working with [ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl/) types with SymPy.jl and Symbolics.jl symbols as coefficients.
- Generation of C-code for filtering with LTI systems.

This package exports the names `s,z` of type `SymPy.Sym` for the Laplace and Z-transform variables. These can be used to build symbolic transfer functions. To build symbolic transfer functions with Symbolics.jl symbols, create `s` or `z` using `using Symbolics; @variables s`.

## Installation
Users typically want to install both ControlSystemsBase and SymbolicControlSystems. ControlSystemsBase contains the basic control-systems functionality, like system types etc., that is used when working with SymbolicControlSystems.
```julia
using Pkg
Pkg.add(["ControlSystemsBase", "SymbolicControlSystems"])
```

## Usage examples
```julia
julia> using ControlSystemsBase, SymbolicControlSystems

julia> @vars w T d # Define (SymPy) symbolic variables
(w, T, d)

julia> h = 0.01; # Sample time

julia> G = tf([w^2], [1, 2*d*w, w^2]) * tf(1, [T, 1])
TransferFunction{Continuous, SisoRational{Sym}}
                      w^2
-----------------------------------------------
T*s^3 + 2*T*d*w + 1*s^2 + T*w^2 + 2*d*w*s + w^2

Continuous-time transfer function model

julia> Gd = tustin(G, h); # Discretize

julia> Sym(G) # Convert a TransferFunction to symbolic expression
                        2                      
                       w                       
───────────────────────────────────────────────
   3    2                   ⎛   2        ⎞    2
T⋅s  + s ⋅(2⋅T⋅d⋅w + 1) + s⋅⎝T⋅w  + 2⋅d⋅w⎠ + w 

julia> ex = w^2 / (s^2 + 2*d*w*s + w^2) # Define symbolic expression
         2       
        w        
─────────────────
           2    2
2⋅d⋅s⋅w + s  + w 

julia> tf(ex) # Convert symbolic expression to TransferFunction
TransferFunction{Continuous, SisoRational{Sym}}
         w^2
---------------------
1*s^2 + 2*d*w*s + w^2

Continuous-time transfer function model

julia> # Replace symbols with numbers
       T_, d_, w_ = 0.03, 0.2, 2.0 # Define system parameters
(0.03, 0.2, 2.0)

julia> Gd_  = sym2num(Gd, h, Pair.((T, d, w), (T_, d_, w_))...)
TransferFunction{Discrete{Float64}, SisoRational{Float64}}
1.4227382019434605e-5*z^3 + 4.2682146058303814e-5*z^2 + 4.2682146058303814e-5*z + 1.4227382019434605e-5
-------------------------------------------------------------------------------------------------------
              1.0*z^3 - 2.705920013658287*z^2 + 2.414628594192383*z - 0.7085947614779404

Sample Time: 0.01 (seconds)
Discrete-time transfer function model
```
Get a Latex-string
```julia
julia> latextf(G)
"\$\\dfrac{1.0w^2}{0.003s^3 + s^2(0.006dw + 1.0) + s(2.0dw + 0.003w^2) + 1.0w^2}\$"
```


## Code generation
The function `code = SymbolicControlSystems.ccode(G::LTISystem)` returns a string with C-code for filtering of a signal through the linear system `G`. All symbolic variables present in `G` will be expected as inputs to the generated function. The transfer-function state is handled by the C concept of `static` variables, i.e., a variable that remembers it's value since the last function invocation. The signature of the generated function `transfer_function` expects all input arguments in alphabetical order, except for the input `u` which always comes first.

Code generation for systems with multiple inputs and outputs (MIMO) is only handled for statespace systems, call `ss(G)` to convert a transfer function to a statespace system.


A usage example follows
```julia
using ControlSystems, SymbolicControlSystems

@vars w T d      # Define symbolic variables
h        = 0.01  # Sample time
G        = tf([w^2], [1, 2*d*w, w^2]) * tf(1, [T, 1])
Gd       = tustin(G, h) # Discretize
code     = SymbolicControlSystems.ccode(Gd, cse=true)
path     = mktempdir()
filename = joinpath(path, "code.c")
outname  = joinpath(path, "test.so")
write(joinpath(path, filename), code)
run(`gcc $filename -lm -shared -o $outname`)

## Test that the C-code generates the same output as lsim in Julia

function c_lsim(u, T, d, w)
    Libc.Libdl.dlopen(outname) do lib
        fn = Libc.Libdl.dlsym(lib, :transfer_function)
        map(u) do u
            @ccall $(fn)(u::Float64, T::Float64, d::Float64, w::Float64)::Float64
        end
    end
end

u    = randn(1,100); # Random input signal 
T_, d_, w_ = 0.03, 0.2, 2.0 # Define system parameters
y    = c_lsim( u,  T_,  d_,  w_); # Filter u through the C-function filter
Gd_  = sym2num(Gd, h, Pair.((T, d, w), (T_, d_, w_))...) # Replace symbols with numeric constants
y_,_ = lsim(ss(Gd_), u); # Filter using Julia
@test norm(y-y_)/norm(y_) < 1e-10
plot([u; y; y_]', lab=["u" "y c-code" "y julia"]) |> display
```

**NOTE:** Numerical accuracy
> The usual caveats for transfer-function filtering applies. High-order transfer functions might cause numerical problems. Consider either filtering through many smaller transfer function in series, or convert the system into a well-balanced statespace system and generate code for this instead. See [lecture notes](http://www.control.lth.se/fileadmin/control/Education/EngineeringProgram/FRTN01/lectures/L11_slides6.pdf) slide 45 and onwards as well as the [ControlSystems docs on numerical accuracy.](https://juliacontrol.github.io/ControlSystems.jl/latest/man/numerical/#Performance-considerations). The function `ControlSystems.ss` converts a transfer function to a statespace system and performs automatic balancing. 


### C-code for gain scheduled systems
The following example writes C-code that interpolates between two linear systems.
The interpolation vector `t` defines the interpolation points.

The system in the example is a double-mass-spring damper, where the inertia of the load is allowed to vary:

```julia
function double_mass_model(Jl) # Inertia load
    Jm = 1  # Inertia motor
    k = 100 # Spring constant
    c0 = 1  # Dampings
    c1 = 1
    c2 = 1
    A = [
        0.0 1 0 0
        -k/Jm -(c1 + c0)/Jm k/Jm c1/Jm
        0 0 0 1
        k/Jl c1/Jl -k/Jl -(c1 + c2)/Jl
    ]
    B = [0, 1/Jm, 0, 0]
    C = [1 0 0 0]
    ss(A,B,C,0)
end

t = [1, 5] # The different inertias in the interpolation
sys = [c2d(double_mass_model(inertia), 0.01) for inertia in t]
SymbolicControlSystems.print_c_array(stdout, sys, t, "mass_spring_damper")
```

This will print C-code functions for the interpolation of each of the system matrices. See the docstring for `print_c_array` for more customization options.