#!/usr/bin/env julia

using EAGOParametricInterval

# write your own tests here
println("Testing Checks and Utilities...")
t = @elapsed include("ParamChk_Tests.jl")
println("done (took $t seconds).")

println("Testing Sparse Preconditioner...")
t = @elapsed include("SparseCntr_Tests.jl")
println("done (took $t seconds).")

println("Testing Parametric Contractors...")
t = @elapsed include("ParametricContractor_Tests.jl")
println("done (took $t seconds).")
