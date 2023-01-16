
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
    τ::Float64 = 5.0      # activity time constant - milliseconds
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


function velocity_input(ωᵢ::OneForm, v::Vector, on_mfld_x::Vector, J::Matrix)
    # ωᵢ(on_mfld_x, J*v) 
    ωᵢ(on_mfld_x, v) 
end


function pushforward(ρ::Function, x::Vector)::Matrix
    J = jacobian(ρ, x)


    # perturb `x` if jacobian has nans
    # η = 0.01
    # while any(isnan.(J))
    #     J = jacobian(ρ, x + rand(size(x)) .* η)
    #     η *= 1.25
    # end
    J[isnan.(J)] .= 0
    return J
end

"""
    step!(simulation::Simulation, x::Vector, v::Vector) 

Step the simulation dynamics given that the "particle" is at `x`
and moving with velocity vector `v`.
"""
function step!(
    simulation::Simulation,
    decoded_x::Vector,
    on_mfld_x::Vector,
    v::Vector;
    s₀ = nothing,
)
    # prep variables
    can = simulation.can
    b₀ = simulation.b₀
    S, W = simulation.S, simulation.W
    Ṡ = simulation.Ṡ .* 0.0
    d = size(S, 2)

    # get velocity input
    J = pushforward(can.C.ρ, decoded_x)
    V = can.α .* map(i -> velocity_input(can.Ω[i], v, on_mfld_x, J), 1:length(can.Ω)) |> vec  # inputs vector of size 2d

    # r(x) = round(x; digits=4)
    # println("\n\n" * string(r.(v)))
    # println(string(r.(V)))

    # update each population with each population's input
    for i = 1:d, j = 1:d
        # get baseline and noise inputs
        input = simulation.η > 0 ? (rand(Float64, size(S, 1)) .* simulation.η) .+ b₀ : b₀

        # enforce initial condition
        isnothing(s₀) || (S[:, i] .*= s₀)

        # get activation
        Ṡ[:, i] .+= W[j] * S[:, j] .+ V[i] .+ input
    end


    # update activity
    simulation.S += (can.σ.(Ṡ) - S) / (simulation.τ)

    # remove bad entries
    droptol!(simulation.S, 0.001)
    droptol!(simulation.Ṡ, 0.001)
    simulation.S = sparse(simulation.S)
    simulation.Ṡ = sparse(simulation.Ṡ)

    return ∑ⱼ(S)  # return the sum of all activations
end




# ---------------------------------------------------------------------------- #
#                                      RUN                                     #
# ---------------------------------------------------------------------------- #

include("history.jl")

function run_simulation(
    simulation::Simulation;
    savefolder::Union{String,Nothing} = nothing,
    savename::String = simulation.can.name * "_sim",
    frame_every_n::Union{Nothing,Int} = 20,   # how frequently to save an animation frame
    fps = 20,
    discard_first_ms = 0,
    s₀ = nothing,
    φ::Union{Function,Nothing} = nothing,
    kwargs...,
)
    savefolder = isnothing(savefolder) ? savename : savefolder

    # setup animation
    N = size(simulation.trajectory.X, 1)
    T = (N + 1) * simulation.dt
    time = 1:simulation.dt:T |> collect
    framen = 1
    anim = Animation()

    # get history to track data
    history = History(simulation, N; discard_first_ms = discard_first_ms, kwargs...)

    # prep some variables
    X̄ = zeros(size(simulation.trajectory.X)) # store decoded variable
    X̄[1, :] = simulation.trajectory.X[1, :]
    decoder_initialized = false
    decoder = nothing

    # do simulation steps and visualize   
    pbar = ProgressBar()
    Progress.with(pbar) do
        job = addjob!(pbar, description = "Simulation", N = N)
        for i = 1:N
            # get activation for bump initialization
            if i > simulation.trajectory.still
                s₀ = nothing
            end

            # get trajectory data
            v = simulation.trajectory.V[min(i+1, N), :]

            # decode manifold bump position
            s̄ = ∑ⱼ(simulation.S)
            if decoder_initialized
                decoded_x, on_mfld_x = decoder(s̄, simulation.can)
            else
                decoded_x, on_mfld_x = simulation.trajectory.X[i, :], decode_peak_location(s̄, simulation.can)
            end
            X̄[i, :] = decoded_x

            # step simulation
            S̄ = step!(simulation, decoded_x, on_mfld_x, v; s₀ = s₀)

            # initialize decoder if necessary
            if (i >= simulation.trajectory.still + 0) && !decoder_initialized
                # prep decoder
                decoder = Decoder(
                    simulation.trajectory.X[i, :],
                    decode_peak_location(S̄, simulation.can),
                    # 1 / simulation.can.offset_size,
                )
                decoder_initialized = true
            end

            # add data to history
            (time[framen] > discard_first_ms) && add!(history, framen, simulation, v)

            # add frame to animation
            isnothing(frame_every_n) || begin
                (i % frame_every_n == 0 || i == 1) &&
                    (time[framen] > discard_first_ms) &&
                    # (framen > simulation.trajectory.still) &&
                    begin
                        # try
                            plot(simulation, time[framen], framen, x, v, X̄, φ)
                        # catch e
                        #     @warn "cacca" framen e
                        # end
                        frame(anim)
                    end
            end

            framen += 1
            update!(job)
        end
    end

    isnothing(frame_every_n) || begin
        gif(anim, savepath(savefolder, savename, "gif"), fps = fps)
    end

    save_simulation_history(
        history,
        savefolder,
        savename * "_" * simulation.can.name * "_history",
    )
    # save_model(
    #     simulation.can,
    #     savefolder,
    #     savename * "_" * simulation.can.name * "_sim_CAN_model",
    #     :CAN,
    # )
    save_data(
        simulation.trajectory.X,
        savefolder,
        savename * "_" * simulation.can.name * "_sim_trajectory_X",
    )
    # save_data(
    #     simulation.trajectory.V,
    #     savefolder,
    #     savename * "_" * simulation.can.name * "_sim_trajectory_V",
    # )
    save_data(X̄, savefolder, savename * "_" * simulation.can.name * "_sim_decoded_X")

    return history, X̄
end
