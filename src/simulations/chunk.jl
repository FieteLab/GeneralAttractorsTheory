abstract type AbstractChunk end


# ---------------------------------------------------------------------------- #
#                                constant chunk                                #
# ---------------------------------------------------------------------------- #
struct ConstantChunk <: AbstractChunk
    duration::Int  # duration in ms
    nframes::Int
    v::Vector{Float64}
end


function ConstantChunk(v::Vector, simulation; duration::Int=500)
    nframes = (Int ∘ floor)(duration / simulation.dt)

    A = simulation.can.A
    @assert length(v) == size(A, 2)  "Chunk vector length $(length(v)) incompatible with A's size $(size(A))"
    return ConstantChunk(duration, nframes, v)
end

# ---------------------------------------------------------------------------- #
#                                 random chunks                                #
# ---------------------------------------------------------------------------- #

struct RandomChunk <: AbstractChunk
    duration::Int  # duration in ms
    nframes::Int
    v::Vector{Vector{Float64}}
end

"""
Pseudo random input combining smoothed white noise with
a sine wave
"""
function RandomChunk(simulation; duration::Int=500, smoothing_window::Int=25)
    nframes = (Int ∘ floor)(duration / simulation.dt)


    smoothing_window = (Int ∘ floor)(smoothing_window / simulation.dt)
    θ = 2π .* clamp.(
        moving_average(
            rand(Float64, nframes), smoothing_window
            ), 
        -1, 1
    ) .+ 2π .* sin.((1:nframes) ./ (150 / simulation.dt))
    
    d = simulation.can.d
    v = if d == 1
        cos.(θ) |> collect
    elseif d == 2
        collect.(zip(sin.(θ), cos.(θ)))
    else
        error("Not implemented for d > 2")
    end

    return RandomChunk(duration, nframes, v)
end