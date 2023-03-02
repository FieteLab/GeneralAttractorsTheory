
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
    trajectory::AbstractTrajectory
    S::SparseMatrixCSC               # n_neurons x 2d - state of all neurons
    W::Union{Nothing, Vector{SparseMatrixCSC} }      # all connection weights
    Ṡ::SparseMatrixCSC
    b₀::Float64 = 1.0       # baseline input activity
    η::Float64 = 0.1       # noise scale
    dt::Float64 = 0.5       # simulation step - milliseconds
    τ::Float64 = 5.0      # activity time constant - milliseconds
end

Base.string(sim::Simulation) = "Simulation of $(sim.can)"
Base.print(io::IO, sim::Simulation) = print(io, string(sim))
Base.show(io::IO, ::MIME"text/plain", sim::Simulation) = print(io, string(sim))

function Simulation(can::AbstractCAN, trajectory::AbstractTrajectory; kwargs...)
    n_pops = hasfield(typeof(can), :offsets) ?  length(can.offsets) : 1

    # initialize activity matrices
    N = *(can.n...)
    S = spzeros(Float64, N, n_pops)
    Ṡ = spzeros(Float64, N, n_pops)

    # get all connection weights
    W = if n_pops > 1
        droptol!.(sparse.(map(x -> Float64.(x), can.Ws)), 0.01)
    else
        nothing
    end

    return Simulation(can = can, trajectory = trajectory, S = S, Ṡ = Ṡ, W = W; kwargs...)
end


# ---------------------------------------------------------------------------- #
#                                     STEP                                     #
# ---------------------------------------------------------------------------- #
∑ⱼ(x) = sum(x, dims = 2) |> vec


"""
    step!(simulation::Simulation, x::Vector, v::Vector) 

Step the simulation dynamics given that the "particle" is at `x`
and moving with velocity vector `v`.
"""
function step!(
    S_tot,
    simulation::Simulation,
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

    # get baseline and noise inputs
    input = simulation.η > 0 ? (rand(Float64, size(S, 1)) .* simulation.η) .+ b₀ : b₀

    # update each population with each population's input
    for i = 1:d
        # enforce initial condition
        isnothing(s₀) || (S[:, i] .*= s₀)

        # get activation
        v_input = can.α * can.Ω[i](on_mfld_x, v)
        Ṡ[:, i] .+= W[i] * S_tot .+ v_input  .+ input
    end

    # update activity
    simulation.S += (can.σ.(Ṡ) - S) / (simulation.τ)

    # remove bad entries
    droptol!(simulation.S, 0.01)
    droptol!(simulation.Ṡ, 0.01)
    simulation.S = sparse(simulation.S)
    simulation.Ṡ = sparse(simulation.Ṡ)

    return ∑ⱼ(simulation.S)  # return the sum of all activations
end



function step!(
    ::Vector,
    simulation::Simulation,
    ::Vector,
    ::Nothing;
    s₀ = nothing,
)
    # prep variables
    can = simulation.can
    b₀ = simulation.b₀
    S, W = simulation.S, simulation.can.W
    Ṡ = simulation.Ṡ .* 0.0


    # get baseline and noise inputs
    input = simulation.η > 0 ? (rand(Float64, size(S, 1)) .* simulation.η) .+ b₀ : b₀

    # enforce initial condition
    isnothing(s₀) || (S .*= s₀)

    # get activation
    Ṡ = W * S .+ input

    # update activity
    simulation.S += (can.σ.(Ṡ) - S) / (simulation.τ)

    # remove bad entries
    droptol!(simulation.S, 0.001)
    droptol!(simulation.Ṡ, 0.001)
    simulation.S = sparse(simulation.S)
    simulation.Ṡ = sparse(simulation.Ṡ)

    return simulation.S  # return the sum of all activations
end


# ---------------------------------------------------------------------------- #
#                                      RUN                                     #
# ---------------------------------------------------------------------------- #

include("history.jl")

function run_simulation(
    simulation::Simulation;
    discard_first_ms = 0,
    s₀ = nothing,
    callbacks::Dict = Dict(),
    kwargs...,
)
    # setup history
    N = size(simulation.trajectory.X, 1)  # number of steps
    T = (N + 1) * simulation.dt
    time = 1:simulation.dt:T |> collect
    framen = 1

    # get history to track data
    history = History(simulation, N; discard_first_ms = discard_first_ms, kwargs...)

    # prep some variables
    X̄ = zeros(size(simulation.trajectory.X)) # store decoded variable
    X̄[1, :] = simulation.trajectory.X[1, :]
    decoder_initialized = false
    decoder = nothing

    # keep track of which callbacks have been invoked already
    cbs_called = Dict(k => false for k in keys(callbacks))


    function do_step(i)
        # get activation for bump initialization
        if i > simulation.trajectory.still
            s₀ = nothing
        end

        # get trajectory data
        v = isnothing(simulation.trajectory.V) ? nothing : simulation.trajectory.V[min(i+1, N), :]

        # decode manifold bump position
        s̄ = ∑ⱼ(simulation.S)
        if decoder_initialized
            decoded_x, on_mfld_x = decoder(s̄, simulation.can)
        else
            decoded_x, on_mfld_x = simulation.trajectory.X[i, :], decode_peak_location(s̄, simulation.can)
        end
        X̄[i, :] = decoded_x 

        # step simulation
        s̄ = step!(s̄, simulation, on_mfld_x, v; s₀ = s₀)

        # initialize decoder if necessary
        if (i >= simulation.trajectory.still + 0) && !decoder_initialized
            _x = decode_peak_location(s̄, simulation.can)

            offset = simulation.can.C.M == simulation.can.C>N ? 0 : simulation.trajectory.X[i, :] .- _x
            # prep decoder
            decoder = Decoder(
                simulation.trajectory.X[i, :],
                _x;
                decoding_offset = offset,
            )
            decoder_initialized = true
        end

        # add data to history
        (time[framen] > discard_first_ms) && add!(history, framen, simulation, v, 
            simulation.trajectory.X[i, :], decoded_x, on_mfld_x
        )

        # call eventual callback functions
        for (name, (cbtime, cb)) in callbacks
            if (time[framen]) >= cbtime && !cbs_called[name] == true
                @debug "Simulation callback $name called at time $(time[framen]), frame $framen" cb name time[framen] framen
                cb(simulation)
                cbs_called[name] = true
            end
        end

        framen += 1
    end

    # do simulation steps and visualize   
    if N > 250  # only show progress bar when its worthwhile
        pbar = ProgressBar()
        Progress.with(pbar) do
            job = addjob!(pbar, description = "Simulation", N = N)
            for i = 1:N
                do_step(i)
                update!(job)
            end
        end
    else
        for i = 1:N
            do_step(i)
        end
    end

    return history, X̄
end
