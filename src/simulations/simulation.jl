
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

function Simulation(can::AbstractCAN; kwargs...)
    # initialize activity matrices
    N = *(can.n...)
    S = spzeros(Float64, N, 2can.d)
    Ṡ = spzeros(Float64, N, 2can.d)

    # get all connection weights
    # W = cat(can.Ws..., dims=3)  # n×n×2d
    W = sparse.(map(x -> Float64.(x), can.Ws))
    droptol!.(W, 0.001)

    return Simulation(can = can, S = S, Ṡ = Ṡ, W = W; kwargs...)
end


# ---------------------------------------------------------------------------- #
#                                     STEP                                     #
# ---------------------------------------------------------------------------- #
∑ⱼ(x) = sum(x, dims = 2) |> vec

function step!(simulation::Simulation, v::Vector{Float64})
    can = simulation.can
    b₀ = simulation.b₀
    A, S, W = simulation.can.A, simulation.S, simulation.W
    Ṡ = simulation.Ṡ

    # get effect of recurrent connectivity & external input
    d = 2simulation.can.d
    B = b₀ .+ A * v  # inputs vector of size 2d
    S̄ = ∑ⱼ(S)  # get the sum of all current activations

    if simulation.η > 0
        η = rand(Float64, size(S, 1), d) .* simulation.η  # get noise input
        for i = 1:d
            Ṡ[:, i] .= W[i] * S̄ .+ B[i] .+ η[i]
        end
    else
        for i = 1:d
            Ṡ[:, i] .= W[i] * S̄ .+ B[i]
        end
    end

    # remove bad entries
    droptol!(simulation.S, 0.001)
    droptol!(simulation.Ṡ, 0.001)

    # update activity
    simulation.S += (can.σ.(Ṡ) - S) / (simulation.τ)
end


# ---------------------------------------------------------------------------- #
#                                    HISTORY                                   #
# ---------------------------------------------------------------------------- #    

""" store simulation data """
mutable struct History
    S::Array                        # activation at each timestep, n_neurons × 2d × T
    Ŝ::Array                        # for averaging
    v::Array                        # input vector at each timestep
    v̂::Array                        # for averaging
    average_over::Int               # average every N frames
    discard::Int                    # discard first N frames
    metadata::Dict                  # can info, sim params, timestamp...
    entry_n::Int                    # keep track of were we should be updating stuff
end


function History(
    simulation::Simulation,
    nframes::Int;
    average_over_ms::Int = 10,
    discard_first_ms = 250,
)
    # see how many frames are we skipping at start
    n_discard = (Int ∘ round)(discard_first_ms / simulation.dt)

    # see over how many frames we average
    average_over = (Int ∘ round)(average_over_ms / simulation.dt)
    keep_frames = (Int ∘ floor)((nframes - n_discard) / average_over)
    keep_frames < 1 && error("Keep frames < 0, reduce discard or increase duration")

    @info "Creating history arrays" size(simulation.S) size(simulation.can.A) keep_frames average_over

    S = Array{Float64}(undef, (size(simulation.S)..., keep_frames))
    Ŝ = Array{Float64}(undef, (size(simulation.S)..., average_over))
    v = Array{Float64}(undef, (size(simulation.can.A, 2), keep_frames))
    v̂ = Array{Float64}(undef, (size(simulation.can.A, 2), average_over))
    @info "Done" size(S) size(Ŝ) size(v) size(v̂)
    metadata = Dict{Symbol,Any}(
        :can => simulation.can.name,
        :n => simulation.can.n,
        :kernel => (string ∘ typeof)(simulation.can.kernel),
        :σ => simulation.can.σ,
        :A => simulation.can.A,
        :b₀ => simulation.b₀,
        :η => simulation.η,
        :dt => simulation.dt,
        :τ => simulation.τ,
        :average_over_ms => average_over_ms,
    )

    @info "Simulation history saving: $(size(S)[end]) frames" "($(round((nframes-n_discard)*simulation.dt; digits=3)) ms tot , averaging every $average_over_ms ms)" "Discarding first $n_discard frames ($discard_first_ms ms)"
    return History(S, Ŝ, v, v̂, average_over, n_discard, metadata, 1)
end

function add!(history::History, framen::Int, simulation::Simulation, v::Vector{Float64})
    framen < history.discard && return  # skip first N frames

    # make sure it fits in history
    framen > size(history.S, 3) * history.average_over + history.discard && begin
        # max_d = size(history.S, 3)*history.average_over+history.discard
        # @warn "Frame too large during add! to history" framen size(history.S) max_d
        return
    end

    # add to "averaging buffer"
    k = size(history.Ŝ, 3)
    Ŝ, v̂ = history.Ŝ, history.v̂

    F = mod(framen, k)
    Ŝ[:, :, F+1] = simulation.S
    v̂[:, F+1] = v

    # update main registry
    if F == 0 || framen == 1
        history.entry_n > size(history.S, 3) && begin
            @warn "Can't append to history, ran out of space"
            return
        end
        history.S[:, :, history.entry_n] = mean(Ŝ; dims = 3)
        history.v[:, history.entry_n] = mean(v̂; dims = 2)
        history.entry_n += 1
    end
end


Base.string(sim::History) = """
    History
    ----------
    can: '$(sim.metadata[:can])'
    nframes: $(size(sim.S)[end])
"""

Base.print(io::IO, sim::History) = print(io, string(sim))
Base.show(io::IO, ::MIME"text/plain", sim::History) = print(io, string(sim))




# ---------------------------------------------------------------------------- #
#                                      RUN                                     #
# ---------------------------------------------------------------------------- #

function run_simulation(
    simulation::Simulation,
    chunks::Vector;
    savename::String = simulation.can.name * "_sim",
    frame_every_n::Union{Nothing,Int} = 20,   # how frequently to save an animation frame
    kwargs...,
)
    @assert eltype(chunks) <: AbstractChunk

    # setup animation
    T = sum(getfield.(chunks, :duration))
    time = 1:simulation.dt:T |> collect
    framen = 1
    anim = Animation()

    # get history to track data
    tot_frames = sum(getfield.(chunks, :nframes))
    history = History(simulation, tot_frames; kwargs...)

    # do simulation steps and visualize
    pbar = ProgressBar()
    Progress.with(pbar) do
        job = addjob!(pbar, description = "Simulation", N = length(time) + 1)
        for chunk in chunks
            for i = 1:chunk.nframes
                v = eltype(chunk.v) == Float64 ? chunk.v : chunk.v[i]
                step!(simulation, v)
                # framen > 50 && break

                # add data to history
                add!(history, framen, simulation, v)

                # add frame to animation
                isnothing(frame_every_n) || begin
                    i % frame_every_n == 0 &&
                        framen < length(time) &&
                        begin
                            plot(simulation, time[framen], v)
                            frame(anim)
                        end
                end

                framen += 1
                update!(job)
            end
        end
    end

    isnothing(frame_every_n) || begin
        @info "saving animation"
        gif(anim, savepath(savename, savename, "gif"), fps = 20)
    end
    save_simulation_history(history, savename, savename)
    return history
end
