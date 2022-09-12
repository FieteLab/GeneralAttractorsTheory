# ---------------------------------------------------------------------------- #
#                                  SIMULATION                                  #
# ---------------------------------------------------------------------------- #

@with_kw_noshow mutable struct Simulation
    can::AbstractCAN
    S::Matrix{Float64}   # n_neurons x 2d - state of all neurons
    Ṡ::Vector{Float64}   # n_neurons - holds activity change for one copy of neurons
    A::Matrix{Float64}   # variable space Φ → ℝᵈ projection matrix
    b₀::Float64 = 0.5    # baseline input activity
    η::Float64  = 0.1    # noise scale
    dt::Float64 = 5      # simulation step - milliseconds
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


# ---------------------------------------------------------------------------- #
#                                     STEP                                     #
# ---------------------------------------------------------------------------- #
function step!(
    simulation::Simulation, v::Vector{Float64}
)

    can = simulation.can
    Ṡ   = simulation.Ṡ
    b₀  = simulation.b₀

    for i in 1:2simulation.can.d
        Ṡ *= 0.0

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
        Ṡ .+= W * sum(simulation.S, dims=2)

        # get effect of bias and velocity input
        Aᵢ = can.offset_directions[i]
        Ṡ .+= b₀ + Aᵢ ⋅ v

        # update S
        simulation.S[:, i] += (can.σ.(Ṡ) - Sᵢ) / simulation.τ
    end
end



# ---------------------------------------------------------------------------- #
#                                      RUN                                     #
# ---------------------------------------------------------------------------- #

function run_simulation(simulation::Simulation, chunks::Vector{SimulationChunk})
        
    # setup animation
    T = sum(getfield.(chunks, :duration))
    time = 1:simulation.dt:T |> collect
    framen = 1
    anim = Animation()

    # do simulation steps and visualize
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

run_simulation(can::AbstractCAN, chunks::Vector{SimulationChunk}) = run_simulation(Simulation(can), chunks)