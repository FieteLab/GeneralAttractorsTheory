
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
Base.show(io::IO, ::MIME"text/plain", sim::Simulation) = print(io, string(sim))



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
    S::Array                # activation at each timestep, n_neurons × 2d × T
    Ŝ::Array                # for averaging
    v::Array                # input vector at each timestep
    v̂::Array                # for averaging
    W::Array{Float64}       # n_neurons × n_neurons × 2d - all connection weights
    average_over::Int       # average every N frames
    discard::Int            # discard first N frames
    metadata::Dict          # can info, sim params, timestamp...
    entry_n::Int            # keep track of were we should be updating stuff
end


function History(simulation::Simulation, nframes::Int; average_over_ms::Int=10, discard_first_ms=250)
    # see how many frames are we skipping at start
    n_discard  = (Int ∘ round)(discard_first_ms / simulation.dt)

    # see over how many frames we average
    average_over = (Int ∘ round)(average_over_ms / simulation.dt)
    keep_frames = (Int ∘ floor)((nframes-n_discard) / average_over)

    S = Array{Float64}(undef, (size(simulation.S)..., keep_frames))
    Ŝ = Array{Float64}(undef, (size(simulation.S)..., average_over))
    v = Array{Float64}(undef, (size(simulation.can.A, simulation.can.d), keep_frames))
    v̂ = Array{Float64}(undef, (size(simulation.can.A, simulation.can.d), average_over))
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
        :average_over_ms=>average_over_ms,
    )

    @info "Simulation history saving: $(size(S)[end]) frames" "($(round((nframes-n_discard)*simulation.dt; digits=3)) ms tot , averaging every $average_over_ms ms)" "Discarding first $n_discard frames ($discard_first_ms ms)"
    return History(S, Ŝ, v, v̂, simulation.W, average_over, n_discard,  metadata, 1)
end

function add!(history::History, framen::Int, simulation::Simulation, v::Vector{Float64})
    framen < history.discard && return  # skip first N frames
    
    # make sure it fits in history
    framen > size(history.S, 3)*history.average_over+history.discard && begin
        @info "skipping" size(history.S, 3) history.average_over framen k*history.average_over+history.discard
        error()
        return
end

    # add to "averaging buffer"
    k = size(history.Ŝ, 3)
    Ŝ, v̂ = history.Ŝ, history.v̂
    
    F = mod(framen, k)
    Ŝ[:, :, F+1] = simulation.S
    v̂[:, F+1] = v

    # update main registry
    if F == 0 || framen==1
        history.entry_n > size(history.S, 3) && begin
            @warn "Can't append to history, ran out of space"
            return
        end
        history.S[:, :, history.entry_n] = mean(Ŝ; dims=3)
        history.v[:, history.entry_n] = mean(v̂; dims=2)
        history.entry_n +=1
    end
end


Base.string(sim::History) = """
    History
    ----------
    can: '$(sim.metadata[:can])'
    nframes: $(size(sim.S)[end])
"""

Base.print(io::IO, sim::History) = print(io, string(sim))
Base.show(io::IO, ::MIME"text/plain", sim::History) = print(io, string(sim))




# ---------------------------------------------------------------------------- #
#                                      RUN                                     #
# ---------------------------------------------------------------------------- #

function run_simulation(
        simulation::Simulation, 
        chunks::Vector;
        savename::String=simulation.can.name*"_sim",
        frame_every_n::Union{Nothing, Int} = 20,   # how frequently to save an animation frame
        kwargs...
    )
    @assert eltype(chunks) <: AbstractChunk
        
    # setup animation
    T = sum(getfield.(chunks, :duration))
    time = 1:simulation.dt:T |> collect
    framen = 1
    anim = Animation()

    # get history to track data
    tot_frames = sum(getfield.(chunks, :nframes))
    history = History(simulation, tot_frames; kwargs...)

    # do simulation steps and visualize
    pbar = ProgressBar()
    Progress.with(pbar) do
        job = addjob!(pbar, description="Simulation",  N=length(time)+1)
        for chunk in chunks
            for i in 1:chunk.nframes
                v = eltype(chunk.v) == Float64 ? chunk.v : chunk.v[i] 
                step!(simulation, v)
                # framen > 300 && break

                # add data to history
                add!(history, framen, simulation, v)

                # add frame to animation
                isnothing(frame_every_n) || begin
                    i % frame_every_n == 0 && framen < length(time) && begin
                        plot(simulation, time[framen], v)
                        frame(anim)
                    end
                end

                framen += 1
                update!(job)
            end
        end
    end
    
    isnothing(frame_every_n) || gif(anim, "test.gif", fps=20)
    save_simulation_history(history, savename)
    return history
end

