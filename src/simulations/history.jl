
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
    v::Array                        # input vector at each timestep
    v̂::Array                        # for averaging
    average_over::Int               # average every N frames
    discard::Int                    # discard first N frames
    metadata::Dict                  # can info, sim params, timestamp...
    entry_n::Int                    # keep track of were we should be updating stuff
    Δt::Number                      # time interval between saved frames. Depends on average over and dt
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
    average_over = average_over_ms > 0 ?  (Int ∘ round)(average_over_ms / simulation.dt) : 1
    keep_frames = (Int ∘ round)((nframes - n_discard) / average_over)
    keep_frames < 1 && error("Keep frames < 0, reduce discard or increase duration")
    Δt = average_over_ms > 0 ? average_over_ms : simulation.dt

    @debug "Creating history arrays" size(simulation.S) keep_frames average_over

    S = Array{Float64}(undef, (size(simulation.S)..., keep_frames))
    Ŝ = Array{Float64}(undef, (size(simulation.S)..., average_over))
    v = Array{Float64}(undef, (size(simulation.trajectory.V, 2), keep_frames))
    v̂ = Array{Float64}(undef, (size(simulation.trajectory.V, 2), average_over))
    
    @debug "Done" size(S) size(Ŝ) size(v) size(v̂)
    metadata = Dict{Symbol,Any}(
        :can => simulation.can.name,
        :cover => simulation.can.C,
        :n => simulation.can.n,
        :kernel => (string ∘ typeof)(simulation.can.kernel),
        :σ => simulation.can.σ,
        :b₀ => simulation.b₀,
        :η => simulation.η,
        :dt => simulation.dt,
        :τ => simulation.τ,
        :average_over_ms => average_over_ms,
    )

    @debug "Simulation history saving: $(size(S)[end]) frames" "($(round((nframes-n_discard)*simulation.dt; digits=3)) ms tot , averaging every $average_over_ms ms)" "Discarding first $n_discard frames ($discard_first_ms ms)"
    return History(S, Ŝ, v, v̂, average_over, n_discard, metadata, 1, Δt)
end

function add!(history::History, framen::Int, simulation::Simulation, v::Vector{Float64})
    framen < history.discard && return  # skip first N frames

    # make sure it fits in history
    framen > size(history.S, 3) * history.average_over + history.discard && begin
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
            # @warn "Can't append to history, ran out of space" F framen
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
