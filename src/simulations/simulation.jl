# ---------------------------------------------------------------------------- #
#                                  SIMULATION                                  #
# ---------------------------------------------------------------------------- #

@with_kw_noshow mutable struct Simulation
    can::AbstractCAN
    S::Matrix{Float64}      # n_neurons x 2d - state of all neurons
    W::Array{Float64}       # n_neurons × n_neurons × 2d - all connection weights
    A::Matrix{Float64}      # 2x × dim(Φ) - variable vel vector projection
    Ṡ::Matrix{Float64}
    b₀::Float64 = 1.0       # baseline input activity
    η::Float64  = 0.1       # noise scale
    α::Float64  = 0.10315   # input scaling
    dt::Float64 = 0.5       # simulation step - milliseconds
    τ::Float64  = 10.0      # activity time constant - milliseconds
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
    # initialize activity matrices
    N = *(can.n...)
    S = zeros(Float64, N, 2can.d)
    Ṡ = zeros(Float64, N, 2can.d)

    # get all connection weights
    W = cat(can.Ws..., dims=3)  # n×n×2d

    # get Φ → neural activity projection mtx
    A = Matrix(hcat(can.offset_directions...)')

    return Simulation(can=can, S=S, Ṡ=Ṡ, W=W, A=A; kwargs...)
end


# ---------------------------------------------------------------------------- #
#                                     STEP                                     #
# ---------------------------------------------------------------------------- #
∑ⱼ(x) = sum(x, dims=2)
function step!(
    simulation::Simulation, v::Vector{Float64}
)
    can         = simulation.can
    b₀, α       = simulation.b₀, simulation.α
    A, S, W     = simulation.A, simulation.S, simulation.W
    Ṡ           = simulation.Ṡ


    # initialize Ṡ with noise
    # Ṡ = rand(Float64, size(S)) * simulation.η
    η() = rand(size(S, 1)) * simulation.η

    # get effect of recurrent connectivity & external input
    d = 2simulation.can.d
    B = b₀ .+ α*(A*v)  # inputs vector of size 2d
    for i in 1:d
        @views Ṡ[:, i] = W[:, :, i] * ∑ⱼ(S) .+ B[i] + η()
        # @views Ṡ[:, i] .+= B[i]
    end

    # update activity
    simulation.S += (can.σ.(Ṡ) - S)/(simulation.τ)
end

# function step!(
#     simulation::Simulation, v::Vector{Float64}
# )

#     can = simulation.can
#     Ṡ   = simulation.Ṡ
#     b₀  = simulation.b₀
#     S_tot = sum(simulation.S, dims=2)

#     for i in 1:2simulation.can.d
#         Ṡ *= 0.0

#         # get random excitatory input 
#         Ṡ += rand(Float64, length(Ṡ)) * simulation.η

#         # get state and connection matrix.
#         Sᵢ::Vector{Float64} = simulation.S[:, i] 
#         W::Matrix{Float64} = can.Ws[i]

#         # get effect of all copies of neurons lattice       
#         Ṡ .+= W * S_tot

#         # get effect of bias and velocity input
#         Aᵢ = can.offset_directions[i]
#         Ṡ .+= b₀ + simulation.α .* (Aᵢ⋅v)

#         # update S
#         simulation.S[:, i] += (can.σ.(Ṡ) - Sᵢ) / simulation.τ
#     end
# end



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
        job = addjob!(pbar, description="Simulation",  N=length(time)+10)
        for chunk in chunks
            for i in 1:chunk.nframes
                step!(simulation, chunk.v)
                framen > 500 && break
                # add frame to animation
                i % 20 == 0 && framen < length(time) && begin
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