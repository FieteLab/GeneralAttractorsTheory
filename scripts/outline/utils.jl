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
function constant_traj_sim(can::SingleCAN, duration, dt, still, τ, b₀)
    # initialize trajectory and simulation
    nframes = (Int ∘ round)(duration / dt)
    trajectory = ConstantTrajectory(
        can;
        T = nframes,
        still = still,
    )

    return Simulation(can, trajectory; η = 0.0, b₀ = b₀, τ = τ)
end

"""
    simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀)

Get and run a simulation for a SingleCAN wwith random initial condition
and constant input trajectory.
"""
function simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀; x₀=nothing)
    x₀, activate = random_init(can; x₀=x₀)
    sim = constant_traj_sim(can, duration, dt, still, τ, b₀)
    h, X = run_simulation(    
        sim;
        discard_first_ms = still,
        average_over_ms = 1,
        s₀ = 1.0 .* activate,
    );
    return h, X
end