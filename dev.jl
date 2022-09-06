using GeneralAttractors

import Base.Iterators: product as ×  # cartesian product
import Distances: Metric, PeriodicEuclidean, pairwise
import StaticArrays: SVector, SA_F32, SMatrix
using Term.Progress
using Plots

# ring attractor CAN - can I do it?


"""
Create a CAN given :
- `n`: the number of neurons in each dimension  
- `ξ`: the position of each neuron, such that ξ(i) is the position of the i-th
    neuron in the neural lattice.
"""
function builder(n::NTuple{N,Int}, ξ::Function, d::Metric, kernel::Kernel) where N
    # check that ξ has the right form
    nargs = first(methods(ξ)).nargs - 1
    @assert N == nargs  "ξ should accept $(N) arguments, accepts: $(nargs)"

    rtype = Base.return_types(ξ, NTuple{N, Int})[1]
    @assert rtype <: AbstractVector "ξ should return a Vector with neuron coordinates, not $rtype"

    # get the index of every neuron in the lattice | Array of size `n`
    lattice_idxs::AbstractArray{NTuple{N, Int}} = ×(map(_n -> 1:_n, n)...) |> collect

    # get the coordinates of every neurons | Array of size `n`
    X::Vector{SVector} = [ξ(idx...) for idx in lattice_idxs]
    @degbug "X" size(X) typeof(X) eltype(X) 

    # get connectivity matrix by applying the kernel function to the pairwise distance mtx
    D =  pairwise(d, hcat(X...))
    @degbug "got pairwise" size(D)

    W = kernel.(D)
    @degbug "Got connectivity" size(W)

    return nothing
end




# ----------------------------------- ring ----------------------------------- #
@info "doing ring"
n = (100,)  # number of neurons in the ring
ξ_r(i::Int)::SVector = SA_F32[2π/n[1]*i]  # neurons coordinates function
d = PeriodicEuclidean([2π])  # distance function


a = 1.0
α = 0.10315
l = 1.0
λ = 13
β = 3/(λ^2)
γ = 1.05 * β

k(x)::Float64 = a * exp(-γ*abs(x)^2) - exp(-β*abs(x)^2)

builder(n, ξ_r, d, Kernel(k));



# # ----------------------------------- torus ---------------------------------- #
# @info "doing torus"
# n = (64, 64)  # number of neurons in the ring
# ξ_t(i::Int, j::Int)::SVector = SA_F32[Float64(i), Float64(j)]  # neurons coordinates function
# d = PeriodicEuclidean([2π, 2π])  # distance function

# @time builder(n, ξ_t, d)