struct SimulationChunk
    duration::Int  # duration in ms
    nframes::Int
    v::Vector
end


function SimulationChunk(v::Vector, simulation; duration::Int=500)
    nframes = (Int ∘ floor)(duration / simulation.dt)

    A = simulation.can.A
    @assert length(v) == size(A, 2)  "Chunk vector length $(length(v)) incompatible with A's size $(size(A))"
    return SimulationChunk(duration, nframes, v)
end

