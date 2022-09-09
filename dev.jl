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
    b₀::Float64 = 0.5    # baseline input activity
    η::Float64  = 0.1    # noise scale
    dt::Float64 = 5    # simulation step - milliseconds
    τ::Float64  = 10.0   # activity time constant - milliseconds
end

Base.string(sim::Simulation) = """
    Simulation
    ----------
    $(sim.can)

    S: $(size(sim.S))
    A: $(size(sim.A))
    η: $(round(sim.η, digits=2))
    dt: $(sim.dt)
    τ: $(sim.τ)
"""

Base.print(io::IO, sim::Simulation) = print(io, string(sim))
Base.show(io::IO, ::MIME"plain/text", sim::Simulation) = print(io, string(sim))




function Simulation(can::AbstractCAN; kwargs...)
    # initialize activity matrix
    N = *(can.n...)
    S = zeros(Float64, N, 2can.d)

    # allocate Ṡ to speed up computation
    Ṡ = zeros(Float64, N)

    # define Φ → can mapping matrix
    A = hcat(can.offset_directions...)
    return Simulation(can=can, S=S, Ṡ=Ṡ, A=A; kwargs...)
end

function Plots.plot(simulation::Simulation, timems, v::Vector; kwargs...)
    can = simulation.can
    plt = plot(;
        title="elapsed: $(round(timems)) ms", 
        # clims=(0.0, 0.05),
        aspect_ratio=:equal, 
        grid=false,
    )

    
    for i in 1:can.d * 2
        S = simulation.S[:, i]
        S = reshape(S, can.n)'

        # get offset position for plotting
        ai = can.offset_directions[i]
        offset = ceil.(ai) .* can.n
        x = collect(1:can.n[2]) .+ offset[1]
        y = collect(1:can.n[1]) .+ offset[2]

        heatmap!(x, y, S)
        x0, y0 = can.n[1]/2, can.n[2]/2
        plot!([x0, x0+v[1]*20], [y0, y0+v[2]*20], lw=6, color=:black, label=nothing)
        scatter!([x0], [y0], ms=8, color=:black, label=nothing)

        # add text with Aᵢ
        Ai = string(round.(ai; digits=3))
        vi = ai ⋅ v
        annotate!(
            [offset[1]+5, x0], 
            [offset[2]+20, y0+20],      
            [
                text(" I: $i\nAi: $Ai\nVi: $(round(vi; digits=4))", :white, :left, 7),
                text("v⃗: $(round.(v; digits=2))", :black, :center, 8),
            ]
        )
    end
    plt
end

# ---------------------------------------------------------------------------- #
#                                     utils                                    #
# ---------------------------------------------------------------------------- #
reset(simulation::Simulation) = begin
    simulation.Ṡ .*= 0.0
end
reset(Ṡ::Vector) = Ṡ * 0.0



# ---------------------------------------------------------------------------- #
#                                     STEP                                     #
# ---------------------------------------------------------------------------- #
function step!(
            simulation::Simulation, v::Vector{Float64}
        )

    can = simulation.can
    Ṡ = simulation.Ṡ
    b₀, A = simulation.b₀, simulation.A

    for i in 1:2simulation.can.d
        # rest Ṡ to 0.0
        reset(Ṡ)

        # get random excitatory input 
        Ṡ += rand(Float64, length(Ṡ)) * simulation.η

        # get state and connection matrix.
        Sᵢ::Vector{Float64} = simulation.S[:, i] 
        W::Matrix{Float64} = can.Ws[i]

        # get effect of all copies of neurons
        # for j in 1:2can.d
        #     j == i || continue
        #     Sⱼ = simulation.S[:, j]
        #     Ṡ .+= W * Sⱼ  # ! ALLOCATIONS HEAVY
        # end        
        weights = ones(4) * .25
        weights[i] = 1.0
        Ṡ .+= (W) * sum(simulation.S * weights, dims=2)

        # get effect of bias and velocity input
        Aᵢ = can.offset_directions[i]
        Ṡ .+= b₀ + Aᵢ ⋅ v

        # update S
        simulation.S[:, i] += (can.σ.(Ṡ) - Sᵢ) / simulation.τ
    end
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


function run(can::AbstractCAN)
    simulation = Simulation(can)
    println(simulation)

    # chunks = [
    #     SimulationChunk(simulation, [cos(0), sin(0)]),
    #     SimulationChunk(simulation, [cos(π/5), sin(π/5)]),
    #     SimulationChunk(simulation, [cos(π/2-π/5), sin(π/2-π/5)]),
    # ]

    chunks = map(
        θ -> SimulationChunk(simulation, [cos(θ), sin(θ)]; duration=250),
        [0, π/4, π/2, 3/4*π, π, -3/4*π, -π/2, -π/4]
    ) |> collect

    # setup animation
    T = sum(getfield.(chunks, :duration))
    time = 1:simulation.dt:T |> collect
    framen = 1
    anim = Animation()

    pbar = ProgressBar()
    Progress.with(pbar) do
        job = addjob!(pbar, description="Simulation",  N=length(time)+1)
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

run(torus_attractor)
# TODO improve performance
# TODO check validity
# TODO add visualization

# TODO test with non-square A
# TODO test with 1D networks