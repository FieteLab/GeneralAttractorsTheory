"""
Estimate attractor dynamics by:

1. loading torus simulation data
2. picking a sample plot
3. for n times, perturn the state in a random direction and measure
    dynamics along the manifold and off the manifold. 
4. to get the on/off manifold dynamics, use local PCA to get the tangent
    space and project the perturbation onto it.

----

Since local PCA depends on NN tree which is hard to 
compute in hiigh dimensions, we use isomap to get a lower dimensional
"""

using ManifoldLearning, MultivariateStats  # needed to be able to load iso
using Distances, FileIO, LinearAlgebra

GENERATE_DATA = false

include("settings.jl")
move_to_datadir(supervisor, "mfld_top")

# --------------------------------- load data -------------------------------- #
# load torus data
tag = "manifold_topology_data"
filters = Dict{Symbol, Any}(
    :tag => tag,
    :can => "torus",
)

@isdefined(X) || begin
    X = load_and_concat_activations(; filters...) # [:, 1:100:end]
end
@isdefined(can) || (can = make_standard_single_torus_can())


# ----------------------------- define functions ----------------------------- #
# write a function to get a random n-dimensional vector of magnitude offset size
function get_random_offset(n::Int)
    offset = randn(n)
    offset ./= norm(offset)
    offset .* offset_size
    offset
end


# write a function to apply the offset to the simulation's state
function apply_offset_to_state!(sim, offset)
    sim.S .+= offset
end


"""
3d plot of a vector at a point
"""
plot_vec!(ax, x, v; kwargs...) = plot!(
    ax, 
    [x[1], x[1] + v[1]], 
    [x[2], x[2] + v[2]];
    # [x[3], x[3] + v[3]];
    kwargs...
)

function vec2mat(v::Vector)::Matrix
    n= length(v)
    reshape(v, (n, 1))
end

# --------------------------------- run sims --------------------------------- #
# get the magnitude of the colums of X to set the perturbation size
X_mags = norm.(eachcol(X))
offset_size = .1 .* mean(X_mags)

# run a simulation applying the stimulust at some point
N = *(can.n...)
x₀ = [0, 3.14]
duration, still = N+250, 25
offset_time = 200
tstart, tend = offset_time-100-still, offset_time+100-still  # time window around stim onset

N_sims = 100

GENERATE_DATA && begin
    tspace_offset_norm = [] # store norm of offset vector in tangent space
    offmld_offset_norm = [] # store norm of offset vector in orthogonal space
    tspace_distances = []  # store distance from pre-stimulus state in tangent space
    offmld_distances = []  # store distance from pre-stimulus state in orthogonal space
end

for simn in 1:N_sims
    GENERATE_DATA || break
    offset = get_random_offset(N)

    @info "running sim $simn/$N_sims"
    h, _ = simulate_constant_traj_random_init(can, 
                duration, dt, still, τ, b₀; 
                x₀=x₀, 
                callbacks=Dict(:offset => (offset_time, 
                    sim -> apply_offset_to_state!(sim, offset)
                    )
                )
    )

    s = h.S[:, 1, 2:end]


    # get the closest points from torus data to the pre-stim trajectory & fit local pca
    pw_d_s = pairwise(
                Euclidean(), 
                vec2mat(s[:, tstart-1]), 
                X
        )[1, :]  # get the closest points based on the time just before the offset stimulus

    distance_theshold = sort(pw_d_s)[501]
    closest = X[:, pw_d_s .< distance_theshold]

    pca = fit(PCA, closest, maxoutdim=400, pratio=1.0)

    # get the tangent space
    tspace = eigvecs(pca)[:, 1:2]

    # project to a space orthogonal to the tangent space (N × N)
    pcs_orthonormal = Matrix(qr(tspace).Q)
    projection_matrix = Matrix{Float64}(I, N, N) - pcs_orthonormal * pcs_orthonormal'

    # @info "Fitted local PCA" s closest pca eigvecs(pca) offset tspace projection_matrix norm(tspace' * offset) norm(projection_matrix * offset)
    
    # store offset vector magnitude
    push!(
        tspace_offset_norm, 
        norm(tspace' * offset)
    )
    push!(
        offmld_offset_norm, 
        norm(projection_matrix * offset)
    )

    # get distance from pre-stimulus state in/out tangent space
    tspace_s = tspace' * s
    tgt = reshape(tspace_s[:, tstart-1], (2, 1))
    push!(tspace_distances, 
        pairwise(Euclidean(), tgt, tspace_s)[1, tstart:tend]
    )

    # project data to orthogonal space and get the dynamics
    off_mfld_s = projection_matrix * s
    tgt = reshape(off_mfld_s[:, tstart-1], (size(off_mfld_s, 1), 1))
    push!(
        offmld_distances,
        pairwise(Euclidean(), tgt, off_mfld_s)[1, tstart:tend]
    )
end


# ---------------------------------- viusals --------------------------------- #

p1 = histogram(tspace_offset_norm, label="tangent space offset vector norm", grid=false)
histogram!(offmld_offset_norm, label="orthogonal space offset vector norm", grid=false)

p2 = plot(
    label = "tangent space distance from pre-stimulus state",
    xlabel = "time (ms)", ylabel = "distance (a.u.)",
    title = "On manifold drift",
    grid=false
)
p3 = plot(
    label = "orthogonal space distance from pre-stimulus state",
    xlabel = "time (ms)", ylabel = "distance (a.u.)",
    title = "Off manifold drift",
    grid=falsegrid=false,
)

for i in 1:N_sims
    plot!(p2, tspace_distances[i], label=nothing, lw=.5, color=:grey, alpha=.75)
    plot!(p3, offmld_distances[i], label=nothing, lw=.5, color=:grey, alpha=.75)
end

# plot mean and std for the distances
tspace_μ = mean(hcat(tspace_distances...), dims=2)
tspace_σ = std(hcat(tspace_distances...), dims=2)
offmld_μ = mean(hcat(offmld_distances...), dims=2)
offmld_σ = std(hcat(offmld_distances...), dims=2)

plot!(p2, tspace_μ, ribbon=tspace_σ, label="mean ± std", color=salmon_dark, lw=3, ribbon_alpha=.5)
plot!(p3, offmld_μ, ribbon=offmld_σ, label="mean ± std", color=indigo, lw=3, ribbon_alpha=.5)


fig = plot(p1, p2, p3, layout=(3, 1), size=(800, 800))


save_plot(supervisor, fig, "on_off_mfld_drift")
fig