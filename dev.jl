using GeneralAttractors
using LinearAlgebra: I, ⋅  # import identity matrix and dot product
# using Tullio
using BenchmarkTools: @benchmark

import GeneralAttractors.Can: AbstractCAN
using Parameters: @with_kw_noshow
using Plots
using Term.Progress


@with_kw_noshow mutable struct Simulation
    can::AbstractCAN
    S::Matrix{Float64}   # n_neurons x 2d - state of all neurons
    Ṡ::Vector{Float64}   # n_neurons - holds activity change for one copy of neurons
    A::Matrix{Float64}   # variable space Φ → ℝᵈ projection matrix
    b₀::Float64 = 0.0    # bias?
    η::Float64  = 0.1    # noise scale
    dt::Float64 = 10.0    # simulation step - milliseconds
    τ::Float64  = 10.0   # activity time constant - milliseconds
end

function Simulation(can::AbstractCAN; A=nothing, kwargs...)
    # initialize activity matrix
    N = *(can.n...)
    S = rand(N, 2can.d)

    # allocate Ṡ to speed up computation
    Ṡ = zeros(N)

    # define Φ → can mapping matrix
    if isnothing(A)
        K = 4  # dimensionality of variable space
        A = Matrix(1.0I, K, can.d)
    else
        K = size(A, 1)
    end

    return Simulation(can=can, S=S, Ṡ=Ṡ, A=A; kwargs...)
end

reset(simulation::Simulation) = begin
    simulation.Ṡ .*= 0.0
end

function step!(
            simulation::Simulation, v::Vector{Float64}
        )

    can = simulation.can
    Ṡ = simulation.Ṡ
    b₀, A = simulation.b₀, simulation.A

    for i in 1:2simulation.can.d
        # get random excitatory input 
        Ṡ += rand(Float64, length(Ṡ)) * simulation.η

        # get state and connection matrix.
        Sᵢ::Vector{Float64} = simulation.S[:, i]  # use view to reduce allocations

        # use view to reduce allocations
        W::Matrix{Float64} = can.Ws[i]

        # get effect of all copies of neurons
        # for j in 1:2can.d
        #     Sⱼ = S[:, j]
        #     Ṡ .+= W * Sⱼ  # ! ALLOCATIONS HEAVY
        # end        
        Ṡ .+= W * sum(simulation.S, dims=2)

        # get effect of bias and velocity input
        a_idx = (Int ∘ ceil)(i / 2)
        Aᵢ::Vector{Float64} = view(A, :, a_idx)
        Ṡ .+= b₀ + Aᵢ ⋅ v  # bias + velocity input

        # update S
        dSdt = (can.σ.(Ṡ) - Sᵢ) / simulation.τ
        # i == 1 && (
        #     @info "up" dSdt[1] Ṡ[1] Sᵢ[1]
        # )
        simulation.S[:, i] += dSdt

        # rest Ṡ
        reset(simulation)
    end
    # println(simulation.S[1:5, 1])

end


function Plots.plot(simulation::Simulation, timems, v::Vector; kwargs...)
    can = simulation.can
    plt = plot(;
        title="elapsed: $(round(timems)) ms", 
        clims=(0.0, 0.1),
        aspect_ratio=:equal, 
    )

    
    for i in 1:can.d * 2
        S = simulation.S[:, i]
        S = reshape(S, can.n)

        offset = can.offset_directions[i] .* can.n
        x = collect(1:can.n[2]) .+ offset[1]
        y = collect(1:can.n[1]) .+ offset[2]

        heatmap!(x, y, S)
        x0, y0 = can.n[1]/2, can.n[2]/2
        plot!([x0, x0+v[1]*20], [y0, y0+v[2]*20], lw=6, color=:black, label=nothing)
        scatter!([x0], [y0], ms=8, color=:black, label=nothing)
    end
    plt
end


# ------------------------------------ run ----------------------------------- #

struct SimulationChunk
    duration::Int  # duration in ms
    nframes::Int
    v::Vector
end


function SimulationChunk(simulation::Simulation, v::Vector; duration::Int=500)
    nframes = (Int ∘ floor)(duration / simulation.dt)
    return SimulationChunk(duration, nframes, v)
end


function run()
    # define Φ → can mapping matrix
    K = 2  # dimensionality of variable space
    A = Matrix(1.0I, K, torus_attractor.d) .* .5

    simulation = Simulation(torus_attractor; A=A)

    chunks = [
        SimulationChunk(simulation, [cos(0), sin(0)]),
        SimulationChunk(simulation, [cos(π/5), sin(π/5)]),
        SimulationChunk(simulation, [cos(π/2-π/5), sin(π/2-π/5)]),
        SimulationChunk(simulation, [0.0, 0.0]; duration=1000),
        SimulationChunk(simulation, [1.0, 0.0]; duration=1000),
        SimulationChunk(simulation, [-1.0, 0.0]; duration=1000),
        SimulationChunk(simulation, [0.0, 1.0]; duration=1000),
        SimulationChunk(simulation, [0.0, -1.0]; duration=1000),
    ]

    # setup animation
    T = sum(getfield.(chunks, :duration))
    time = 1:simulation.dt:T |> collect
    framen = 1
    anim = Animation()

    pbar = ProgressBar()
    Progress.with(pbar) do
        job = addjob!(pbar, description="Simulation",  N=length(time))
        for chunk in chunks
            for i in 1:chunk.nframes
                step!(simulation, chunk.v)

                # add frame to animation
                i % 1 == 0 && begin
                    plot(simulation, time[framen], chunk.v)
                    frame(anim)
                end

                framen += 1
                update!(job)
            end
        end
    end

    gif(anim, "test.gif", fps=20)
end

run()
# TODO improve performance
# TODO check validity
# TODO add visualization