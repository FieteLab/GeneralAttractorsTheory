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

GENERATE_DATA = true

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

# project data to 400d dimensional space
d400_pca = @isdefined(d400_pca) ? d400_pca : fit(PCA, X, maxoutdim=4096, pratio=1.0)
d400_P = projection(d400_pca)[:, 1:400]'
d400_x = @isdefined(d400_x) ? d400_x : d400_P * X
@isdefined(pw_d_x) || (pw_d_x = pairwise(Euclidean(), d400_x))

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

N_sims = 3

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


    # project trajectory and offset to 400d space  
    d400_s = (d400_P * h.S[:, 1, :])[:, 2:end]

    # get the closest points from torus data to the pre-stim trajectory & fit local pca
    pw_d_s = pairwise(
                Euclidean(), 
                Matrix(reshape(d400_s[:, tstart-1], (400, 1))), 
                d400_x
        )[1, :]  # get the closest points based on the time just before the offset stimulus

    distance_theshold = sort(pw_d_s)[501]
    closest = d400_x[:, pw_d_s .< distance_theshold]

    pca = fit(PCA, closest, maxoutdim=400, pratio=1.0)
    @info "aa" d400_s closest pca eigvecs(pca)
     
    @info "fulld_pca" fulld_pca
    tspace = Matrix(fulld_pca[:, 1:2]')
    ortho_space = Matrix(fulld_pca[:, 3:end]')

    @assert size(tspace) == (2, N) "tspace is wrong size: $(size(tspace)) instead of $((2, N))"
    @assert size(ortho_space) == (N-2, N) "ortho_space is wrong size: $(size(ortho_space)) instead of $((N-2, N))"
 
    # get norm of offset vector in/out tangent space
    push!(tspace_offset_norm, norm(tspace * offset))
    push!(offmld_offset_norm, norm(ortho_space * offset))

    # get distance from pre-stimulus state in/out tangent space
    tspace_s = tspace * d400_s
    tgt = reshape(tspace_s[:, tstart-1], (2, 1))
    push!(tspace_distances, 
        pairwise(Euclidean(), tgt, tspace_s)[1, tstart:tend]
    )

    # project data to orthogonal space and get the dynamics
    off_mfld_s = ortho_space * d400_s
    tgt = reshape(off_mfld_s[:, 1], (398, 1))
    push!(
        offmld_distances,
        pairwise(Euclidean(), tgt, off_mfld_s)[1, tstart:tend]
    )
end


# ---------------------------------- viusals --------------------------------- #

p1 = histogram(tspace_offset_norm, label="tangent space offset vector norm")
histogram!(offmld_offset_norm, label="orthogonal space offset vector norm")

p2 = plot(
    label = "tangent space distance from pre-stimulus state",
    xlabel = "time (ms)", ylabel = "distance (a.u.)",
    title = "On manifold drift",
)
p3 = plot(
    label = "orthogonal space distance from pre-stimulus state",
    xlabel = "time (ms)", ylabel = "distance (a.u.)",
    title = "Off manifold drift",
)

for i in 1:N_sims
    plot!(p2, tspace_distances[i], label=nothing, lw=.5, color=:grey)
    plot!(p3, offmld_distances[i], label=nothing, lw=.5, color=:grey)
end

fig = plot(p1, p2, p3, layout=(3, 1), size=(800, 800))


