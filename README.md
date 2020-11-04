# SymbolicControlSystems

[![Build Status](https://github.com/baggepinnen/SymbolicControlSystems.jl/workflows/CI/badge.svg)](https://github.com/baggepinnen/SymbolicControlSystems.jl/actions)
[![Coverage](https://codecov.io/gh/baggepinnen/SymbolicControlSystems.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/baggepinnen/SymbolicControlSystems.jl)




## Code generation
The function `code = SymbolicControlSystems.ccode(G::LTISystem)` returns a string with C-code for filtering of a signal through the linear system `G`. A usage example follows
```julia
using ControlSystems, SymbolicControlSystems

@vars w T d # Define symbolic variables
h        = 0.01
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

u    = randn(100); # Random input signal 
T_, d_, w_ = 0.03, 0.2, 2.0 # Define system parameters
y    = c_lsim( u,  T_,  d_,  w_); # Filter u through the C-function filter
Gd_  = sym2num(Gd, h, Pair.((T, d, w), (T_, d_, w_))...) # Replace symbols with numeric constants
y_,_ = lsim(Gd_, u); # Filter using Julia
@test norm(y-y_)/norm(y_) < 1e-10
plot([u y y_], lab=["u" "y c-code" "y julia"]) |> display
```