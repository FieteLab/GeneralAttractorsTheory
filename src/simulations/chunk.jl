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
    v::Vector{Vector}
end


function RandomChunk(simulation; duration::Int=500)
    nframes = (Int ∘ floor)(duration / simulation.dt)

    d = size(simulation.can.A, 2)

    rnd = 2π .* [
        [clamp(x-0.5, -1, 1)] 
        for x in moving_average(
            rand(Float64, seqlen), 11
            )
            ]

    v = zip(sin.(rnd), cos.(rnd)) |> collect
    return RandomChunk(duration, nframes, v)
end