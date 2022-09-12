struct SimulationChunk
    duration::Int  # duration in ms
    nframes::Int
    v::Vector
end


function SimulationChunk(v::Vector; dt::Float64=0.5, duration::Int=500)
    nframes = (Int ∘ floor)(duration / dt)
    return SimulationChunk(duration, nframes, v)
end

