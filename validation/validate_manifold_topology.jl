"""
Validate `GlobalAttractors.ManifoldAnalysis` by analyzing 
data from known manifolds
"""

using Parameters
using GeneralAttractors
using GeneralAttractors.Analysis
using Plots
using Term.Progress
using Statistics

# ---------------------------------------------------------------------------- #
#                                MANIFOLD TYPES                                #
# ---------------------------------------------------------------------------- #

abstract type AbstractPointManifold end

@with_kw struct Torus <: AbstractPointManifold
    d::Int           = 2
    size::Vector     = [2π, 2π]
    ϕ::Function      = (θ₁, θ₂; R=0.75, r=0.25) -> [
        (R + r * cos(θ₁)) * cos(θ₂),
        (R + r * cos(θ₁)) * sin(θ₂),
        r * sin(θ₁),
    ]
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
    N::Int=5000,
    η::Float64=.1
)::Matrix{Float64}

    # sample points on the manifold domain
    X = rand(m.d, N) .* m.size

    # embed in ℝ³
    ϕ(x) =  m.ϕ(x...)   
    M = hcat(ϕ.(eachcol(X))...)
    M .+= η .* rand(size(M)...) .- η/2 # the - is required to center the noise offset
end


# ---------------------------------------------------------------------------- #
#                                   ANALYSIS                                   #
# ---------------------------------------------------------------------------- #
params = AnalysisParameters(
    tda_threshold           = 3,
    tda_downsample_factor   = 1,
    tda_dim_max             = 2,
)
manifold = Torus()
M =  generate_manifold_pointcloud(manifold)
@info "Created manifold" size(M)

# info save an animation
# animate_3d_scatter(M, "test", "M")

# ------------------------------------ tda ----------------------------------- #
# tda = estimate_manifold_topology(M, params)
# TODO make sure we can extract the topology from the barcode


# ------------------------- intrinsic dimensionality ------------------------- #
function test_intrinsic_dimensionality(manifold::AbstractPointManifold)
    N, F = collect(0:.025:1), collect(0.5:.1:5)
    D = []

    pbar = ProgressBar()
    job = addjob!(pbar, N=length(N)*length(F))
    Progress.with(pbar) do
        for η in N, f in F
            params = AnalysisParameters(
                intrinsic_d_data_fraction_threshold=f,
            )
            M = generate_manifold_pointcloud(manifold; η=η)

            @info typeof(M) typeof(params)
            estimate_intrinsic_dimensionality(M, params; verbose=false)
            return
            push!(D, estimate_intrinsic_dimensionality(M, params; verbose=false) |> mean)

            update!(pbar)
        end
    end

    contourf(
        N, F, reshape(D, (length(N), length(F))),
        xlabel="Noise η",
        ylabel="frac x ∈ Uᵢ",
        title = "Estimated intrinsic dimensionality"
    )
end

test_intrinsic_dimensionality(manifold)
