
# ---------------------------------------------------------------------------- #
#                                MANIFOLD TYPES                                #
# ---------------------------------------------------------------------------- #

abstract type AbstractPointManifold end

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
        p1 
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


# -------------------------------- higher dim -------------------------------- #
abstract type HigherDimensionalManifold <: AbstractPointManifold end

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

function generate_manifold_pointcloud(m::HigherDimensionalManifold; N=1000, η=.1)
    M = typeof(m)(m.k; N=N, η=η).pts
    @assert size(M, 2) == N size(M)
    @assert size(M, 1) == m.k size(M)
    
    M
end


# ----------------------------------- utils ---------------------------------- #
function project_to_ℝ³(M::Matrix)
    size(M, 1) <= 3 && return M

    pca = fit(PCA, M; maxoutdim=3)
    return predict(pca, M)

end

