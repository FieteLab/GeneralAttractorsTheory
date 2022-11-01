import ToeplitzMatrices: Circulant
import SparseArrays: SparseMatrixCSC, sparse

# ------------------------------ linear algebra ------------------------------ #
""" transpose function """
function ᵀ end

ᵀ(M::Matrix)::Matrix = M' |> collect
ᵀ(M::SparseMatrixCSC)::SparseMatrixCSC = M' |> collect |> sparse
ᵀ(M::Circulant)::Circulant = M' |> Circulant


# ----------------------------------- misc ----------------------------------- #
"""
    lerp(i::Int, n::Int, x₀::Float64, x₁::Float64)::Float64

Linear interpolation between two values (xᵢ) 
given an index `i ∈ [1, n]`.
"""
function lerp(i::Int, n::Int, x₀, x₁)::Float64
    p = (i - 1) / (n - 1)
    x₀ * (1 - p) + x₁ * p
end


# --------------------------------- smoothing -------------------------------- #
"""
    moving_average(A::AbstractArray, m::Int)

Moving average smoothing.
"""
function moving_average(A::AbstractArray, m::Int)
    out = similar(A)
    R = CartesianIndices(A)
    Ifirst, Ilast = first(R), last(R)
    I1 = m ÷ 2 * oneunit(Ifirst)
    for I in R
        n, s = 0, zero(eltype(out))
        for J = max(Ifirst, I - I1):min(Ilast, I + I1)
            s += A[J]
            n += 1
        end
        out[I] = s / n
    end
    return out
end



# ----------------------------------- sizes ---------------------------------- #
"""
    bounding_box(m::Matrix)::NamedTuple

Min and Max along each dimension (row) of a matrix with points cloud data.
Rows are dimensions, columns are points.
"""
function bounding_box(m::Matrix)::NamedTuple
    (; min = minimum.(eachrow(m)), max = maximum.(eachrow(m)))
end

"""
    bounding_box_size(m::Matrix)::Float64

Magnitude of the diagonal of the bounding box of a points cloud
"""
function bounding_box_size(m::Matrix)::Float64
    bb = bounding_box(m)
    norm(bb.max - bb.min)
end
