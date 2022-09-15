using GeneralAttractors
using GeneralAttractors.Analysis
using Plots
using Statistics
using NearestNeighbors
using Term
using MultivariateStats
using Parameters

sim_name = "torus_sim"  # name of the simulation being anlyzed

# TODO collate manifold analysis functions into a neat pipeline



function generate_manifold_pointcloud(
    m::AbstractPointManifold; 
    N::Int=1000,
    η::Float64=.1
)::Matrix{Float64}

    # sample points on the manifold domainN = 1000
    X = rand(m.d, N) .* m.size

    # embed in ℝ³
    ϕ(x) =  m.ϕ(x...)
    M = hcat(ϕ.(eachcol(X))...)
    M .+= η .* rand(size(M)...) .- η/2 # the - is required to center the noise offset
end

T = Torus() |> generate_manifold_pointcloud
# scatter(T[1, :], T[2, :], T[3, :], ylim=[-1, 1], xlim=[-1, 1], zlim=[-1, 1])

animate_3d_scatter(T)