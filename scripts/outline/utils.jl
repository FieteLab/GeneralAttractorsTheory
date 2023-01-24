import GeneralAttractors.Can: SingleCAN

"""
Get a random initial condition for a CAN
"""
function random_init(can; x₀=nothing)
    x₀ = something(x₀, rand(can.C.N))
    d = map(i -> can.metric(x₀, can.X[:, i]), 1:size(can.X, 2))
    activate = zeros(length(d))
    activate[d.<0.5] .= 1
    return x₀, activate
end

"""
    constant_traj_si(can, duration, dt, still, τ, b₀)

Get a sim for a single CAN with a constant trajectory.
"""
function constant_traj_sim(can::SingleCAN, duration, dt, still, τ, b₀; η=0.1)
    # initialize trajectory and simulation
    nframes = (Int ∘ round)(duration / dt)
    trajectory = ConstantTrajectory(
        can;
        T = nframes,
        still = still,
    )

    return Simulation(can, trajectory; b₀ = b₀, τ = τ, η=0.1)
end

"""
    simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀)

Get and run a simulation for a SingleCAN wwith random initial condition
and constant input trajectory.
"""
function simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀; x₀=nothing, 
    η=0.0,
    )
    x₀, activate = random_init(can; x₀=x₀)
    sim = constant_traj_sim(can, duration, dt, still, τ, b₀; η=η)
    h, X = run_simulation(    
        sim;
        discard_first_ms = still,
        average_over_ms = 1,
        s₀ = 1.0 .* activate,
    );
    return h, X
end



# --------------------------------- topology --------------------------------- #

"""
Load a subset of the data from the supervisor with multiple simulations and concatenate
activations. 
"""
function load_and_concat_activations(; filters...)
    _, data = ProjectSupervisor.fetch(supervisor; filters...)
    
    # stack activations over time
    X = hcat(
        map(
            d -> d["S"][:, 1, end-5:end], data
        )...
    )
    @info "Loaded $(length(data)) simulations." X
    return X
end


"""
do PCA dimensionality reduction
"""
do_pca(X, params) = Analysis.pca_dimensionality_reduction(X, params)[2]

"""
do Isomap dimensionality reduction
"""
do_isomap(X, params) = isomap_dimensionality_reduction(X, params)[2]
