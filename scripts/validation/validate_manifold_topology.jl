"""
Validate `GlobalAttractors.ManifoldAnalysis` by analyzing 
data from known manifolds
"""

using Parameters
using Plots
using Term.Progress
using Statistics: mean
using LinearAlgebra: norm
using LinearAlgebra: ⋅
using MultivariateStats
using NearestNeighbors
import Manifolds: Sphere as MSphere
import Manifolds: uniform_distribution

using GeneralAttractors
using GeneralAttractors.Analysis
import GeneralAttractors: savepath
import GeneralAttractors: bounding_box_size, bounding_box

include("manifolds.jl")

function validate_manifold(manifold::AbstractPointManifold)
    @info "Validating" manifold
    # generate ranges for noise and local PCA bounding_box_size
    # get the size of the bounding box
    M = generate_manifold_pointcloud(manifold)
    σ = bounding_box_size(M)

    # get noise values
    noise = range(0, σ, length=5)./σ |> collect

    # local PCA neighborhood size
    local_pca_nsize = collect(range(2, 500, length=20))
    local_pca_npoints = (Int ∘ round).(collect(range(500, 5000, length=25)))

    # iterate over noise levels and perform: plot, TDA, local PCA
    plots = []
    for (n, η) in enumerate(noise)
        # plot manifold points
        @info "Doing step $n/$(length(noise))"
        M = generate_manifold_pointcloud(manifold; η=η, N=1000) |> project_to_ℝ³
        push!(plots, scatter(eachrow(M)..., label=nothing, color="black", ms=2, 
            xlim=[-1.5, 1.5],
            ylim=[-1.5, 1.5],
            zlim=[-1.5, 1.5],   
            title="Noise scale: $(round(η; digits=2))",
            xlabel="PC1", ylabel="PC2", zlabel="PC3",
            camera=(30, 40), size=(200, 200)
        ))

        # do local PCA
        @info "doing local PCA"
        D = zeros(20, 25)
        for (i, nsize) in enumerate(local_pca_nsize), (j, npoints) in enumerate(local_pca_npoints)
            m = generate_manifold_pointcloud(manifold; η=η, N=npoints) 

            params = AnalysisParameters(
                intrinsic_d_neighborhood_size=nsize,
                intrinsic_d_nseeds=200,
            )
            d = estimate_intrinsic_dimensionality(m, params)
            D[i, j] = mean(d)
        end

        # plot results
        push!(plots, 
            heatmap(
                local_pca_npoints, local_pca_nsize, D,
                xlabel="N manifold points",
                ylabel="neighborhood scale",
                title = n==1 ? "Estimated d" : nothing,
                clims=(manifold.d-1.5, manifold.d+1.5),
                levels = 3,
                lc="black",
                c=:bwr,
                colorbar_title="d error"
            )
        )

        # do TDA
        @info "Doing TDA"
        M = generate_manifold_pointcloud(manifold; η=η, N=2000)
        params = AnalysisParameters(
            tda_threshold           = 3,
            tda_downsample_factor   = 1,
            tda_dim_max             = manifold.d,
        )
        _, tda_plot = estimate_manifold_topology(M, params)
        push!(plots, tda_plot)
    end

    @info "plotting"
    p = plot(
            plots...,
            size=(1400, 1400), 
            layout=(5, 3)
    )
    savefig(
        savepath(
        "manifold_validation", string(manifold), "png"
        )
    )
    return plots
end




# ---------------------------------------------------------------------------- #
mflds = (
    # ProjectedManifold(10, Sphere()),
    ProjectedManifold(10, Plane()),
    ProjectedManifold(10, Ring()),
    ProjectedManifold(10, Torus()),
    ProjectedManifold(10, Cylinder()),
    Torus(),
    Plane(),
    Ring(),
    Cylinder(),
    Sphere()
)
for manifold in mflds
    validate_manifold(manifold)
end

