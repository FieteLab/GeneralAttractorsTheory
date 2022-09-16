"""
Validate `GlobalAttractors.ManifoldAnalysis` by analyzing 
data from known manifolds
"""

using Parameters
using GeneralAttractors
using GeneralAttractors.Analysis
using Plots
using Term.Progress
using Statistics: mean
using LinearAlgebra: norm
using MultivariateStats
using NearestNeighbors

import GeneralAttractors: bounding_box_size, bounding_box

# ---------------------------------------------------------------------------- #
#                                MANIFOLD TYPES                                #
# ---------------------------------------------------------------------------- #

abstract type AbstractPointManifold end

@with_kw struct Plane <: AbstractPointManifold
    d::Int           = 2
    size::Vector     = [1, 1]
    ϕ::Function      = (x, y) -> [
        x, y, x*y
    ]
end

@with_kw struct Torus <: AbstractPointManifold
    d::Int           = 2
    size::Vector     = [2π, 2π]
    ϕ::Function      = (θ₁, θ₂; R=1.0, r=0.75) -> [
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
        sin(p)
    ]
end


# -------------------------------- higher dim -------------------------------- #
abstract type HigherDimensionalManifold <: AbstractPointManifold end

struct Sphere_ℝᵈ <: HigherDimensionalManifold
    d::Int
    k::Int
    pts::Matrix  # d x N

    function Sphere_ℝᵈ(k::Int; N=1000, η=0.0)
        u = rand(-1:.1:1, (k, N*k))

        radius = map(
            x -> abs(norm(x) - 1), eachcol(u)
        )
        u = u[:, radius .< η + 0.01]

        # u = mapslices(
        #     x -> x/norm(x), u; dims=1
        # )

        # get noise
        # η != 0.0 && begin
        #     η = η/2
        #     u .+= η .* rand(size(u)...) .- η/2 
        # end
        new(2, k, u)
    end
end






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
    X = rand(m.d, N*3) .* m.size

    # embed in ℝ³
    ϕ(x) =  m.ϕ(x...)   
    M = hcat(ϕ.(eachcol(X))...)
    η = η/2
    M .+= η .* rand(size(M)...) .- η/2 # the - is required to center the noise offset
end

function generate_manifold_pointcloud(m::HigherDimensionalManifold; N=1000, η=.1)
    typeof(m)(m.k; N=N*m.k, η=η).pts
end


# ----------------------------------- utils ---------------------------------- #
function project_to_ℝ³(M::Matrix)
    size(M, 1) <= 3 && return M

    pca = fit(PCA, M; maxoutdim=3)
    return predict(pca, M)

end


# ---------------------------------------------------------------------------- #
#                                   VISUALIZATION                              #
# ---------------------------------------------------------------------------- #
function plot_manifold_vs_noise(manifold::AbstractPointManifold)
    # get the size of the bounding box
    M = generate_manifold_pointcloud(manifold)
    bb = bounding_box(M)
    σ = bounding_box_size(M)

    # get noise values
    noise = range(0, σ/3, length=9)./σ |> collect

    # plot manifolds
    plt = []
    @track for η in noise 
        M = generate_manifold_pointcloud(manifold; η=η, N=500) |> project_to_ℝ³
        
        p = plot()
        scatter!(eachrow(M)..., label=nothing, color="black", ms=2)
        plot!(
            title="Noise: $(round(η*100; digits=2))%",
            # xlim=[-1, 1], ylim=[-1, 1], zlim=[-1, 1],
            camera=(30, 40)
        )

        push!(plt, p)
    end
    plot(plt...; size=(1000, 800)) |> display
end


"""
Visualize points in the neighborhood of a given point 
selected using NN
"""
function visualize_nn(manifold::AbstractPointManifold; K=70)

    M = generate_manifold_pointcloud(manifold; N=2000)
    nntree = KDTree(M)

    seed = M[:, rand(1:size(M, 2))]

    # U = inrange(nntree, seed, ϵ)
    U, _ = knn(nntree, seed, K) 
    M̄ = project_to_ℝ³(M)

    p = scatter(eachrow(M̄)..., color=:black, alpha=.8, ms=5,
        xlim=[-1, 1], ylim=[-1, 1], zlim=[-1, 1], camera=(90, 0),
        label=nothing,
    )
    scatter!(eachrow(M̄[:, U])..., color=:red, ms=5, label="in radius")

    plot!(title="K:$K, n points in nn $(length(U))")
    display(p)
end

# ---------------------------------------------------------------------------- #
#                                   ANALYSIS                                   #
# ---------------------------------------------------------------------------- #
# params = AnalysisParameters(
#     tda_threshold           = 3,
#     tda_downsample_factor   = 1,
#     tda_dim_max             = 2,
# )
# manifold = Torus()
# M =  generate_manifold_pointcloud(manifold)
# @info "Created manifold" size(M)

# info save an animation
# animate_3d_scatter(M, "test", "M")



# ------------------------------------ tda ----------------------------------- #
# tda = estimate_manifold_topology(M, params)
# TODO make sure we can extract the topology from the barcode


# ------------------------- intrinsic dimensionality ------------------------- #

function nanmean(x) 
    x̄ = filter(!isnan,x)

    return if length(x̄) > 0
        mean(x̄)
    else
        NaN
    end
end

function test_intrinsic_dimensionality(manifold::AbstractPointManifold; K=3000)
    # get size of bounding box of clean manifold
    M = generate_manifold_pointcloud(manifold)
    σ = bounding_box_size(M)


    n_samples = 10
    # get noise scale (fraction of bb)
    N = collect(range(0, 0.1σ, length=n_samples)./σ)

    # K in knn
    F = collect(range(10, 100, length=n_samples))

    D = []
    pbar = ProgressBar()
    job = addjob!(pbar, N=length(N)*length(F))
    Progress.with(pbar) do
        for η in N, f in F
            params = AnalysisParameters(
                intrinsic_d_bbox_fraction_threshold=f,
                intrinsic_d_npoints=50,
            )

            mfld = generate_manifold_pointcloud(manifold; η=η, N=K)
            d = estimate_intrinsic_dimensionality(mfld, params; verbose=true)
            push!(D, d |> nanmean)
            Progress.update!(job)
        end
    end

    contourf(
        N, F, reshape(D, (length(N), length(F))),
        xlabel="noise scale",
        ylabel="neighborhood scale",
        title = "Estimated intrinsic dimensionality | K=$K",
        clims=(0, 2*manifold.d),
        levels = 4*manifold.d,
        lc="black",
        c=:bwr
    ) |> display
end




# ---------------------------------------------------------------------------- #
#                                      run                                     #
# ---------------------------------------------------------------------------- #
# manifold = Sphere_ℝᵈ(10)
manifold = Cylinder()

# visuals
plot_manifold_vs_noise(manifold)
# visualize_nn(manifold)


# intrinsic dimensionality
# test_intrinsic_dimensionality(manifold)



# ------------------------------------- x ------------------------------------ #
fraction_variance_explained(pca::PCA) = principalvars(pca) ./ var(pca) * 100

M = generate_manifold_pointcloud(manifold; N=10000, η=0.01)
@warn size(M) 


# build nearest neighbor tree
# nntree = BruteTree(M; leafsize=2)
nntree = KDTree(M; leafsize=2)

# sample random points on the manifold
seeds = M[:, rand(1:size(M, 2), 10)]

# get neighborhoods
Us, _ = knn(nntree, seeds, 200) 
# Us = inrange(nntree, seeds, .3)

# for each neighborhood fit PCA and get the number of PCs
p = plot()
for U in Us
    pca_model = fit(PCA, M[:, U]; pratio=.999999);
    plot!(fraction_variance_explained(pca_model), lw=2, label=nothing)
end

display(p)

# TODO work on computing the elbow
# TODO summary plot for each manifold as a fn of noise