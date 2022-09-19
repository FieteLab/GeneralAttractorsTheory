
# ---------------------------------------------------------------------------- #
#                                MANIFOLD TYPES                                #
# ---------------------------------------------------------------------------- #

abstract type AbstractPointManifold end

Base.string(m::AbstractPointManifold) = string(typeof(m))

@with_kw struct Plane <: AbstractPointManifold
    d::Int           = 2
    size::Vector     = [1, 1]
    ϕ::Function      = (x, y) -> [
        2x-1, 2y-1, x*y
    ]
end

@with_kw struct Torus <: AbstractPointManifold
    d::Int           = 2
    size::Vector     = [2π, 2π]
    ϕ::Function      = (θ₁, θ₂; R=2.0, r=1.0) -> [
        (R + r * cos(θ₁)) * cos(θ₂),
        (R + r * cos(θ₁)) * sin(θ₂),
        r * sin(θ₁),
    ]
end


@with_kw struct Sphere <: AbstractPointManifold
    d::Int           = 2
    size::Vector     = [2π, 2π]
    ϕ::Function      = (p0, p1) -> [
        sin(p0) * cos(p1), 
        sin(p0) * sin(p1), 
        cos(p0)
    ]
end

@with_kw struct Cylinder <: AbstractPointManifold
    d::Int           = 2
    size::Vector     = [2π, 1]
    ϕ::Function      = (p0, p1) -> [
        sin(p0),
        cos(p0),
        2p1 - 1 
    ]
end



@with_kw struct Ring <: AbstractPointManifold
    d::Int           = 1
    size::Vector     = [2π]
    ϕ::Function      = (p) -> [
        sin(p),
        cos(p),
        cos(p)
    ]
end


# --------------------------------- projected -------------------------------- #
abstract type ProjectedPointManifold <: AbstractPointManifold end

Base.string(m::ProjectedPointManifold) = "$(string(m.mfld)) embedded in ℝᵏ (k=$(m.k))"

struct ProjectedManifold <: ProjectedPointManifold
    d::Int
    k::Int
    mfld::AbstractPointManifold

    ProjectedManifold(k, mfld::AbstractPointManifold) = new(2, k, mfld)
end


# -------------------------------- higher dim -------------------------------- #
abstract type HigherDimensionalManifold <: AbstractPointManifold end

Base.string(m::HigherDimensionalManifold) = "$(string(typeof(m))) with d=$(m.d)"

struct Sphere_ℝᵈ <: HigherDimensionalManifold
    d::Int
    k::Int
    pts::Matrix  # d x N

    function Sphere_ℝᵈ(k::Int; N=1000, η=0.0)
        sphere = MSphere(k-1)
        u = hcat(rand(uniform_distribution(sphere), N)...)

        u = mapslices(
            x -> x ./ norm(x), u; dims=1
        )
        
        # get noise
        η != 0.0 && begin
            # η = η/2
            u .+= η .* rand(size(u)...) .- η/2 
        end

        new(k-1, k, u)
    end
end

# struct NTorus <: HigherDimensionalManifold
#     d::Int
#     k::Int
#     pts::Matrix  # d x N

#     function Sphere_ℝᵈ(k::Int; N=1000, η=0.0)
#         torus = MTorus(k-1)
#         u = hcat(rand(uniform_distribution(sphere), N)...)
        
#         # get noise
#         η != 0.0 && begin
#             η = η/2
#             u .+= η .* rand(size(u)...) .- η/2 
#         end

#         new(2, k, u)
#     end
# end







# ----------------------------- generate manifold ---------------------------- #
"""
    generate_manifold_pointcloud(
        m::AbstractPointManifold; 
        N::Int=1000,
        η::Float64=.1
    )::Matrix{Float64}

Construct a pointcloud in ℝ³ with data sampled from a
manifold's domain and with added noise using the manifold's
embedding function ϕ. 
"""
function generate_manifold_pointcloud(
    m::AbstractPointManifold; 
    N::Int= 250,
    η::Float64=.1
)::Matrix{Float64}

    # sample points on the manifold domain
    X = rand(m.d, N) .* m.size

    # embed in ℝ³
    ϕ(x) =  m.ϕ(x...)   
    M = hcat(ϕ.(eachcol(X))...)
    η = η/2
    M .+= η .* rand(size(M)...) .- η/2 # the - is required to center the noise offset
end

"""
Get a <4 dimensional manifold and project it onto a 
k-dimensional space
"""
function generate_manifold_pointcloud(m::ProjectedPointManifold; N=1000, η=.1)::Matrix{Float64}
    M = generate_manifold_pointcloud(m.mfld; N=N, η=η)
    P = gram_schmidt(rand(m.k, 3))
    return P * M
end


function generate_manifold_pointcloud(m::HigherDimensionalManifold; N=1000, η=.1)
    M = typeof(m)(m.k; N=N, η=η).pts
    @assert size(M, 2) == N size(M)
    @assert size(M, 1) == m.k size(M)
    
    M
end


# ----------------------------------- utils ---------------------------------- #
function project_to_ℝ³(M::Matrix)
    size(M, 1) <= 3 && return M

    pca = fit(PCA, M; maxoutdim=3, pratio=0.99999999)
    return predict(pca, M)

end

"""
gram schmidth orthogonalization and normalization
of a projection matrix
"""
function gram_schmidt(P::Matrix)::Matrix
    basis = []
    for v in eachrow(P)
        w = isempty(basis) ? v : v .- sum(map(
                                        b -> v ⋅ b, basis
                                    )) 
        push!(basis, w / norm(w))
    end
    return Matrix(hcat(basis...)')
end
