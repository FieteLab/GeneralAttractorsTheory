
# ---------------------------------------------------------------------------- #
#                                  SIMULATION                                  #
# ---------------------------------------------------------------------------- #

"""
    @with_kw_noshow mutable struct Simulation
        can::AbstractCAN
        S::Matrix{Float64}             
        W::Vector{SparseMatrixCSC}      
        Ṡ::Matrix{Float64}
        b₀::Float64 = 1.0      
        η::Float64  = 0.1      
        dt::Float64 = 0.5      
        τ::Float64  = 10.0     
    end

Holds information necessary for running a simulation.

## Arguments
- `can`: CAN network instance
- `S`: placeholder for the state (activation) of all neurons. Shape: n_neurons × 2d
- `W`: holds all network connectivity info as a vector of 2d-many sparse matrices (n_neurons × n_neurons each)
- `Ṡ`: place holder for dS/dt array.
"""
@with_kw_noshow mutable struct Simulation
    can::AbstractCAN
    trajectory::Trajectory
    S::SparseMatrixCSC               # n_neurons x 2d - state of all neurons
    W::Vector{SparseMatrixCSC}       # all connection weights
    Ṡ::SparseMatrixCSC
    b₀::Float64 = 1.0       # baseline input activity
    η::Float64 = 0.1       # noise scale
    dt::Float64 = 0.5       # simulation step - milliseconds
    τ::Float64 = 10.0      # activity time constant - milliseconds
end

Base.string(sim::Simulation) = "Simulation of $(sim.can)"
Base.print(io::IO, sim::Simulation) = print(io, string(sim))
Base.show(io::IO, ::MIME"text/plain", sim::Simulation) = print(io, string(sim))

function Simulation(can::AbstractCAN, trajectory::Trajectory; kwargs...)
    n_pops = length(can.offsets)

    # initialize activity matrices
    N = *(can.n...)
    S = spzeros(Float64, N, n_pops)
    Ṡ = spzeros(Float64, N, n_pops)

    # get all connection weights
    W = sparse.(map(x -> Float64.(x), can.Ws))
    droptol!.(W, 0.001)

    return Simulation(can = can, trajectory = trajectory, S = S, Ṡ = Ṡ, W = W; kwargs...)
end


# ---------------------------------------------------------------------------- #
#                                     STEP                                     #
# ---------------------------------------------------------------------------- #
∑ⱼ(x) = sum(x, dims = 2) |> vec

"""
    step!(simulation::Simulation, x::Vector, v::Vector) 

Step the simulation dynamics given that the "particle" is at `x`
and moving with velocity vector `x`.
"""
function step!(simulation::Simulation, x::Vector, v::Vector; s₀=nothing)
    can = simulation.can
    b₀ = simulation.b₀
    S, W = simulation.S, simulation.W
    Ṡ = simulation.Ṡ

    # get effect of recurrent connectivity & external input
    d = size(S, 2)
    V = vec(
        map(
            # ωᵢ -> ωᵢ(x, v) / norm(ωᵢ(x)), 
            ωᵢ -> ωᵢ(x, v), 
            simulation.can.Ω
            )
    )  # inputs vector of size 2d
    # r(x) = round(x, digits=2)
    # rand() > .75 && println(r.(100 .* v), " "^10, r.(V))

    S̄ = ∑ⱼ(S)  # get the sum of all current activations
    !isnothing(s₀) && (S̄ .*= s₀)

    for i = 1:d
        if simulation.η > 0
            η = rand(Float64, size(S, 1), d) .* simulation.η  # get noise input
            Ṡ[:, i] .= W[i] * S̄ .+ V[i] .+ η[i] .+ b₀
        else
            Ṡ[:, i] .= W[i] * S̄ .+ V[i] .+ b₀
        end
    end

    # remove bad entries
    droptol!(simulation.S, 0.001)
    droptol!(simulation.Ṡ, 0.001)

    # update activity
    simulation.S += (can.σ.(Ṡ) - S) / (simulation.τ)
    return S̄
end




# ---------------------------------------------------------------------------- #
#                                      RUN                                     #
# ---------------------------------------------------------------------------- #

include("history.jl")

function run_simulation(
    simulation::Simulation;
    savename::String = simulation.can.name * "_sim",
    frame_every_n::Union{Nothing,Int} = 20,   # how frequently to save an animation frame
    fps = 20,
    discard_first_ms = 0,
    s₀ = nothing,
    activation_steps=20,
    kwargs...,
)

    # setup animation
    N = size(simulation.trajectory.X, 1)
    T = (N + 1) * simulation.dt
    time = 1:simulation.dt:T |> collect
    framen = 1
    anim = Animation()

    # get history to track data
    history = History(simulation, N; discard_first_ms = discard_first_ms, kwargs...)

    # do simulation steps and visualize
    pbar = ProgressBar()
    X̄ = zeros(size(simulation.trajectory.X))
    X̄[1, :] = simulation.trajectory.X[1, :]
    Progress.with(pbar) do
        job = addjob!(pbar, description = "Simulation", N = N)
        for i = 1:N
            # get activation for bump initialization
            if i > activation_steps
                s₀ = nothing
            end

            # step simulation
            x = simulation.trajectory.X[i, :]
            v = simulation.trajectory.V[i, :]
            # S̄ = step!(simulation, X̄[i, :], v; s₀=s₀)

            x̂ = i < 10 ? simulation.trajectory.X[i, :] : X̄[i-1, :]
            S̄ = step!(simulation, x̂, v; s₀=s₀)

            # decode manifold bump position
            X̄[i, :] = decode_peak_location(S̄, simulation.can)

            # add data to history
            add!(history, framen, simulation, v)

            # add frame to animation
            isnothing(frame_every_n) || begin
                (i % frame_every_n == 0 || i == 1) &&
                    (time[framen] > discard_first_ms) &&
                    begin
                        plot(simulation, time[framen], framen, x, v, X̄)
                        frame(anim)
                    end
            end

            framen += 1
        update!(job)
        end
    end

    isnothing(frame_every_n) || begin
        gif(anim, savepath(savename, savename, "gif"), fps = fps)
    end

    # save_simulation_history(history, savename, savename)
    # save_model(simulation.can, savename, "sim_CAN_model", :CAN)
    # save_data(simulation.trajectory.X, savename, "sim_trajectory_X")
    # save_data(simulation.trajectory.V, savename, "sim_trajectory_V")
    return history, X̄
end
