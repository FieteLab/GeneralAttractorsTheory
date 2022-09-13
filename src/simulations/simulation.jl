
# ---------------------------------------------------------------------------- #
#                                  SIMULATION                                  #
# ---------------------------------------------------------------------------- #

@with_kw_noshow mutable struct Simulation
    can::AbstractCAN
    S::Matrix{Float64}      # n_neurons x 2d - state of all neurons
    W::Array{Float64}       # n_neurons × n_neurons × 2d - all connection weights
    Ṡ::Matrix{Float64}
    b₀::Float64 = 1.0       # baseline input activity
    η::Float64  = 0.1       # noise scale
    dt::Float64 = 0.5       # simulation step - milliseconds
    τ::Float64  = 10.0      # activity time constant - milliseconds
    frame_every_n::Int = 20  # how frequently to save an animation frame
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

    return Simulation(can=can, S=S, Ṡ=Ṡ, W=W; kwargs...)
end


# ---------------------------------------------------------------------------- #
#                                     STEP                                     #
# ---------------------------------------------------------------------------- #
∑ⱼ(x) = sum(x, dims=2) |> vec
function step!(
    simulation::Simulation, v::Vector{Float64}
)   
    can         = simulation.can
    b₀          = simulation.b₀
    A, S, W     = simulation.can.A, simulation.S, simulation.W
    Ṡ           = simulation.Ṡ

    η() = rand(size(S, 1)) * simulation.η

    # get effect of recurrent connectivity & external input
    d = 2simulation.can.d
    B = b₀ .+ A*v  # inputs vector of size 2d

    for i in 1:d
        @views Ṡ[:, i] .= W[:, :, i] * ∑ⱼ(S) .+ B[i] .+ η()
    end

    # update activity
    simulation.S += (can.σ.(Ṡ) - S)/(simulation.τ)
end


# ---------------------------------------------------------------------------- #
#                                    HISTORY                                   #
# ---------------------------------------------------------------------------- #    

""" store simulation data """
mutable struct History
    S::Array        # activation at each timestep, n_neurons × 2d × T
    v::Array        # input vector at each timestep
    W::Array{Float64}       # n_neurons × n_neurons × 2d - all connection weights
    metadata::Dict  # can info, sim params, timestamp...
end


function History(simulation::Simulation, nframes::Int)
    S = Array{Float64}(undef, (size(simulation.S)..., nframes))
    v = Array{Float64}(undef, (size(simulation.can.A, 2), nframes))
    metadata = Dict{Symbol, Any}(
        :can=>simulation.can.name,
        :n=>simulation.can.n,
        :kernel=>(string ∘ typeof)(simulation.can.kernel),
        :σ=>simulation.can.σ,
        :A=>simulation.can.A,
        :b₀=>simulation.b₀,
        :η=>simulation.η,
        :dt=>simulation.dt,
        :τ=>simulation.τ,
    )

    return History(S, v, simulation.W, metadata)
end

function add!(history::History, framen::Int, simulation::Simulation, v::Vector{Float64})
    @assert size(history.v, 2) >= framen "Attempted to instert data for frame $framen and history length $(size(history.v, 2))"
    history.S[:, :, framen] = simulation.S
    history.v[:, framen] = v
end






# ---------------------------------------------------------------------------- #
#                                      RUN                                     #
# ---------------------------------------------------------------------------- #

function run_simulation(
        simulation::Simulation, 
        chunks::Vector;
        savename::String=simulation.can.name*"_sim"
    )
    @assert eltype(chunks) <: AbstractChunk
        
    # setup animation
    T = sum(getfield.(chunks, :duration))
    time = 1:simulation.dt:T |> collect
    framen = 1
    anim = Animation()

    # get history to track data
    tot_frames = sum(getfield.(chunks, :nframes))
    history = History(simulation, tot_frames)

    # do simulation steps and visualize
    pbar = ProgressBar()
    Progress.with(pbar) do
        job = addjob!(pbar, description="Simulation",  N=length(time)+1)
        for chunk in chunks
            for i in 1:chunk.nframes
                v = eltype(chunk.v) isa Number ? chunk.v : chunk.v[i] 
                step!(simulation, v)
                # framen > 300 && break

                # add data to history
                add!(history, framen, simulation, v)

                # add frame to animation
                i % simulation.frame_every_n == 0 && framen < length(time) && begin
                    plot(simulation, time[framen], v)
                    frame(anim)
                end

                framen += 1
                update!(job)
            end
        end
    end
    
    gif(anim, "test.gif", fps=20)
    save_simulation_history(history, savename)
end

