abstract type AbstractChunk end


# ---------------------------------------------------------------------------- #
#                                constant chunk                                #
# ---------------------------------------------------------------------------- #
struct ConstantChunk <: AbstractChunk
    duration::Int  # duration in ms
    nframes::Int
    v::Vector{Float64}
    dt::Float64
end


function ConstantChunk(v::Vector, simulation; duration::Int=500)
    nframes = (Int ∘ floor)(duration / simulation.dt)

    A = simulation.can.A
    @assert length(v) == size(A, 2)  "Chunk vector length $(length(v)) incompatible with A's size $(size(A))"
    return ConstantChunk(duration, nframes, v, simulation.dt)
end

# ---------------------------------------------------------------------------- #
#                                 random chunks                                #
# ---------------------------------------------------------------------------- #

struct RandomChunk <: AbstractChunk
    duration::Int  # duration in ms
    nframes::Int
    v::Vector{Vector{Float64}}
    dt::Float64
end

"""
Pseudo random input combining smoothed white noise with
a sine wave
"""
function RandomChunk(simulation; duration::Int=500, smoothing_window::Int=40, μ₀=2.0, σ=1.0)
    nframes = (Int ∘ floor)(duration / simulation.dt)


    smoothing_window = (Int ∘ floor)(smoothing_window / simulation.dt)

    # create a kinda-random angle vector by summing multiple sine waves and some noise
    η = moving_average(
                rand(Float64, nframes), 25
                )
    θ = sin.((1:nframes) ./ (100 / simulation.dt)) .- 
                sin.((1:nframes) ./ (75 / simulation.dt)) .+ 
                sin.((1:nframes) ./ (125 / simulation.dt)) .+
                sin.((1:nframes) ./ (512 / simulation.dt)) .+
                sin.((1:nframes) ./ (50 / simulation.dt)) .+
                η
    

    μ = μ₀ .+ (moving_average(
        rand(Float64, nframes), smoothing_window
        ) .- 0.5 ) .* σ
    
    d = simulation.can.d
    v = if d == 1
        cos.(θ) .* μ |> collect
    elseif d == 2
        collect.(zip(sin.(θ), cos.(θ))) .* μ
    else
        error("Not implemented for d > 2")
    end

    return RandomChunk(duration, nframes, v, simulation.dt)
end


function Plots.plot(chunk::RandomChunk; kwargs...)
    time = 0:chunk.dt:(chunk.nframes * chunk.dt - chunk.dt)
    d = length(chunk.v[1])
    plt = if d == 1
        y = map(
            x->x[1], chunk.v
        )
        plot(time, y, label=nothing; kwargs...)
    else
        μ = norm.(chunk.v)
        p1 = plot(time, μ, label="magnitude")
        p2 = plot()
        for i in 1:d
            y = map(
                x->x[i], chunk.v
            )
            plot!(time, y, label=nothing; kwargs...)
        end

        θ = acos.(clamp.(
            map(
            x->x[1], chunk.v
        ),-1, 1))
        p3 = plot(θ, label="θ")
        plot(p1, p2, p3)
    end
    display(plt)
end