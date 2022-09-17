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
import Manifolds: Sphere as MSphere
import Manifolds: Torus as MTorus
import Manifolds: uniform_distribution


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
    ϕ::Function      = (θ₁, θ₂; R=1.0, r=0.25) -> [
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
    X = rand(m.d, N*3) .* m.size

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

function project_to_ℝ⁵(M::Matrix)
    size(M, 1) <= 5 && return M

    pca = fit(PCA, M; maxoutdim=5, pratio=0.9999999)
    return predict(pca, M)

end

# ---------------------------------------------------------------------------- #
#                                   VISUALIZATION                              #
# ---------------------------------------------------------------------------- #

function plot!_boundingbox(M::Matrix)

    l, h = bounding_box(M)
    xl, yl, zl = l
    xh, yh, zh = h

    plot!(
        [xl,xh], [yl, yl], [zl, zl], lw=5, color="red", label=nothing,
        xlim=[xl-0.1, xh+0.1],
        ylim=[yl-0.1, yh+0.1],
        zlim=[zl-0.1, zh+0.1],
    )
    plot!(
        [xl,xh], [yh, yh], [zl, zl], lw=5, color="red", label=nothing
    )
    plot!(
        [xl,xh], [yl, yl], [zh, zh], lw=5, color="red", label=nothing
    )
    plot!(
        [xl,xh], [yh, yh], [zh, zh], lw=5, color="red", label=nothing
    )

    plot!(
        [xl,xl], [yl, yh], [zh, zh], lw=5, color="red", label=nothing
    )
    plot!(
        [xh,xh], [yl, yh], [zh, zh], lw=5, color="red", label=nothing
    )
    plot!(
        [xl,xl], [yl, yh], [zl, zl], lw=5, color="red", label=nothing
    )
    plot!(
        [xh,xh], [yl, yh], [zl, zl], lw=5, color="red", label=nothing
    )

    plot!(
        [xl,xl], [yl, yl], [zl, zh], lw=5, color="red", label=nothing
    )
    plot!(
        [xh,xh], [yl, yl], [zl, zh], lw=5, color="red", label=nothing
    )
    plot!(
        [xl,xl], [yh, yh], [zl, zh], lw=5, color="red", label=nothing
    )
    plot!(
        [xh,xh], [yh, yh], [zl, zh], lw=5, color="red", label=nothing
    )
    
end


function plot_manifold_vs_noise(manifold::AbstractPointManifold)
    # get the size of the bounding box
    M = generate_manifold_pointcloud(manifold)
    σ = bounding_box_size(M)

    # get noise values
    noise = range(0, 2σ, length=9)./σ |> collect

    # plot manifolds
    plt = []
    @track for η in noise 
        M = generate_manifold_pointcloud(manifold; η=η, N=500) |> project_to_ℝ³
        
        p = plot()
        plot!_boundingbox(M)
        scatter!(eachrow(M)..., label=nothing, color="black", ms=2)

        bb = bounding_box(M)
        Δ = mean(bb.max .- bb.min)
        plot!(
            title="Noise scale: $(round(η; digits=2))\nBbox scale $(round(Δ; digits=2))",
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
        xlim=[-1.1, 1.1], ylim=[-1.1, 1.1], zlim=[-1.1, 1.1], camera=(90, 0),
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

function test_intrinsic_dimensionality(manifold::AbstractPointManifold; K=10_000)
    # get size of bounding box of clean manifold
    M = generate_manifold_pointcloud(manifold)
    σ = bounding_box_size(M)
    
    # get noise scale (fraction of bb)
    nN, nF = 10, 20
    N = collect(range(0, 1σ, length=nN)./σ)

    # K in knn
    F = collect(range(2, 100, length=nF))

    # estimate intrinsic dimensionality over all combinations
    D = []
    pbar = ProgressBar()
    job = addjob!(pbar, N=length(N)*length(F))
    D = zeros(nN, nF)
    Progress.with(pbar) do
        for (i, η) in enumerate(N)
            m̂ = generate_manifold_pointcloud(manifold; η=η, N=K) # |> project_to_ℝ⁵

            for (j, f) in enumerate(F)
                params = AnalysisParameters(
                    intrinsic_d_bbox_fraction_threshold=f,
                    intrinsic_d_npoints=50,
                )
                d = estimate_intrinsic_dimensionality(m̂, params)
                D[i, j] = mean(d)
                Progress.update!(job)
            end
        end
    end


    contourf(
        F, N, D,
        ylabel="noise scale",
        xlabel="neighborhood scale",
        title = "'$(string(typeof(manifold)))' estimated d | K=$K",
        # clims=(0, 2manifold.d),
        clims=(manifold.d-3, manifold.d+3),
        levels = 4*manifold.d,
        lc="black",
        c=:bwr,
        colorbar_title="estimated d"
    ) |> display

    # heatmap(D) |> display
end




# ---------------------------------------------------------------------------- #
#                                      run                                     #
# ---------------------------------------------------------------------------- #
manifold = Sphere_ℝᵈ(5)
# manifold = Sphere()

# visuals
# plot_manifold_vs_noise(manifold)
# visualize_nn(manifold)


# intrinsic dimensionality
test_intrinsic_dimensionality(manifold)

