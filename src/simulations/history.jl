
# ---------------------------------------------------------------------------- #
#                                    HISTORY                                   #
# ---------------------------------------------------------------------------- #    

""" 
Store simulation data 

In addition to storing data at each frame, history also stores a running
average of all data.
"""
mutable struct History
    S::Array                        # activation at each timestep, n_neurons × 2d × T
    Ŝ::Array                        # for averaging
    v::Union{Nothing, Array}                        # input vector at each timestep
    v̂::Union{Nothing, Array}                        # for averaging
    x_M::Array                       # on mfld position at each timestep
    x̂_M::Array                       # for averagin
    x_M_decoded::Array                       # on mfld position at each timestep
    x̂_M_decoded::Array                       # for averagin
    x_N::Array                       # on mfld position at each timestep
    x̂_N::Array                       # for averagin
    average_over::Int               # average every N frames
    discard::Int                    # discard first N frames
    metadata::Dict                  # can info, sim params, timestamp...
    entry_n::Int                    # keep track of were we should be updating stuff
    Δt::Number                      # time interval between saved frames. Depends on average over and dt
end

function Base.Dict(h::History)
    keys = fieldnames(History)
    vals = map(k -> getfield(h, k), keys)
    return Dict([string(k) => v for (k,v) in zip(keys, vals)]...)
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
    average_over = average_over_ms > 0 ? (Int ∘ round)(average_over_ms / simulation.dt) : 1
    keep_frames = (Int ∘ round)((nframes - n_discard) / average_over)
    keep_frames < 1 && error("Keep frames < 0, reduce discard or increase duration")
    Δt = average_over_ms > 0 ? average_over_ms : simulation.dt

    @debug "Creating history arrays" size(simulation.S) keep_frames average_over nframes

    S = Array{Float64}(undef, (size(simulation.S)..., keep_frames))
    Ŝ = Array{Float64}(undef, (size(simulation.S)..., average_over))
    if isnothing(simulation.trajectory.V)
        v, v̂ = nothing, nothing
    else
        v = Array{Float64}(undef, (size(simulation.trajectory.V, 2), keep_frames))
        v̂ = Array{Float64}(undef, (size(simulation.trajectory.V, 2), average_over))
    end

    d_m, d_n = length(simulation.can.C.M.xmin), length(simulation.can.C.N.xmin)
    x_M = Array{Float64}(undef, (d_m, keep_frames))  # bump position on neural manifold
    x̂_M = Array{Float64}(undef, (d_m, average_over))
    x_M_decoded = Array{Float64}(undef, (d_m, keep_frames))  # bump position on neural manifold
    x̂_M_decoded = Array{Float64}(undef, (d_m, average_over))
    x_N = Array{Float64}(undef, (d_m, keep_frames))  # decoded position on M manifold
    x̂_N = Array{Float64}(undef, (d_m, average_over))
    
    @debug "Done" size(S) 
    metadata = Dict{Symbol,Any}(
        # :can => simulation.can.name,
        # :cover => simulation.can.C,
        # :n => simulation.can.n,
        # :kernel => (string ∘ typeof)(simulation.can.kernel),
        # :σ => simulation.can.σ,
        # :b₀ => simulation.b₀,
        # :η => simulation.η,
        # :dt => simulation.dt,
        # :τ => simulation.τ,
        :average_over_ms => average_over_ms,
    )

    @debug "Simulation history saving: $(size(S)[end]) frames" "($(round((nframes-n_discard)*simulation.dt; digits=3)) ms tot , averaging every $average_over_ms ms)" "Discarding first $n_discard frames ($discard_first_ms ms)"
    return History(S, Ŝ, v, v̂, x_M, x̂_M, x_M_decoded, x̂_M_decoded, x_N, x̂_N, average_over, n_discard, metadata, 1, Δt)
end

function add!(history::History, framen::Int, simulation::Simulation, v::Union{Vector, Nothing}, x_M::Vector, x_M_decoded::Vector, x_N::Vector)
    framen < history.discard && return  # skip first N frames

    # make sure it fits in history
    framen > size(history.S, 3) * history.average_over + history.discard && begin
        return
    end

    # add to "averaging buffer"
    k = size(history.Ŝ, 3)
    Ŝ, v̂, x̂_M, x̂_M_decoded, x̂_N = history.Ŝ, history.v̂, history.x̂_M, history.x̂_M_decoded, history.x̂_N

    F = mod(framen, k)
    Ŝ[:, :, F+1] = simulation.S
    isnothing(v) || (v̂[:, F+1] = v)
    x̂_M[:, F+1] = x_M
    x̂_M_decoded[:, F+1] = x_M_decoded
    x̂_N[:, F+1] = x_N

    # update main registry
    if F == 0 || framen == 1
        history.entry_n > size(history.S, 3) && begin
            # @warn "Can't append to history, ran out of space" F framen
            return
        end
        history.S[:, :, history.entry_n] = mean(Ŝ; dims = 3)
        isnothing(v) || (history.v[:, history.entry_n] = mean(v̂; dims = 2))
        history.x_M[:, history.entry_n] = mean(x̂_M; dims = 2)
        history.x_M_decoded[:, history.entry_n] = mean(x̂_M_decoded; dims = 2)
        history.x_N[:, history.entry_n] = mean(x̂_N; dims = 2)
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
