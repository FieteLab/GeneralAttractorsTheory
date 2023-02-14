import ForwardDiff
using Term.Progress
import FluxTraining: Callback, EpochEnd, Phase, Read, StepEnd, TrainingPhase, EpochBegin, ValidationPhase
using UnicodePlots

# ---------------------------------------------------------------------------- #
#                                     DATA                                     #
# ---------------------------------------------------------------------------- #

"""
    generate_groundtruth_data(trajectory)

Given a trajectory, generate ground-truth data on what `ω`,
the velocity-dependent input to the network should be. 

to test:
    x̄, _ =  test()


    plot(eachcol(trajectory.X)..., lw=4, color=:black)
    plot!(eachrow(hcat(x̄...))..., lw=2, color=:red)

"""
function generate_groundtruth_data(trajectory, warmup)
    S, can_decoder = initialize_can(warmup, trajectory)    
    x̄, Ω = [], []
    for t in 1:size(trajectory.X, 1)
        i = argmax(S)
        ∂w∂x, ∂w∂y = get_W_partial_derivatives(i)

        x = trajectory.X[t, :]
        ẋ = trajectory.V[t, :]
        J = jacobian(x, can)

        v = J*ẋ
        ω = α * v[1] * ∂w∂x + α * v[2] * ∂w∂y
        Ṡ = (can.W * S .+ ω .+ b₀) .|> can.σ
        S += (Ṡ - S)/τ

        push!(x̄, can_decoder(S, can)[1])
        push!(Ω, ω)
    end

    return x̄, Ω
end




# ---------------------------------------------------------------------------- #
#                                      CAN                                     #
# ---------------------------------------------------------------------------- #

"""
jacobian of the cover map of a can
"""
function jacobian(x, can)
    J = ForwardDiff.jacobian(can.C.ρ, x)
    J[isnan.(J)] .= 0
    return J
end


"""
initialize CAN state and decoder at end of warmup
"""
function initialize_can(warmup, trajectory)
    # initialize can state at end of warmup
    S = warmup.S[:, 1, end]

    # initialize S -> X decoder
    x̂ = decode_peak_location(S, can)
    decoder = Decoder(
        trajectory.X[1, :],
        x̂;
        decoding_offset = trajectory.X[1, :] .- x̂,
    )

    return S, decoder
end


o = zeros(can_n, can_n)

"""
    get_W_partial_derivatives

Get the partial derivatives of the weight matrix 
along each direction at neuron `i`.
"""
function get_W_partial_derivatives(i)
    o .*= 0
    Δx = can.X[1, :] .- can.X[1, i]
    wx = can.kernel(Δx)
    ∂w∂x = begin
        o[2:end, :] = diff(reshape(wx, (can_n, can_n)); dims=1)
        vcat(o...)
    end

    o .*= 0
    Δy = can.X[2, :] .- can.X[2, i]
    wy = can.kernel(Δy)
    ∂w∂y = begin
        o[:, 2:end] = diff(reshape(wy, (can_n, can_n)); dims=2)
        vcat(o...)
    end

    return ∂w∂x, ∂w∂y
end


# ---------------------------------------------------------------------------- #
#                                   LEARNING                                   #
# ---------------------------------------------------------------------------- #

# --------------------------- progress bar callback -------------------------- #
mutable struct ProgressCB <: Callback
    progress::ProgressBar
    n_steps::Int
    epochs_job::ProgressJob
    steps_job::Union{Nothing, ProgressJob}
end

function ProgressCB(n_steps::Int, n_epochs::Int)
    pbar = ProgressBar(expand=false, transient=true, columns=:detailed)
    start!(pbar)

    job = addjob!(pbar; N = n_epochs, description="Training run")

    return ProgressCB(pbar, n_steps,job, nothing)
end

function FluxTraining.on(
    event::EpochBegin,
    phase::TrainingPhase,
    cb::ProgressCB,
    learner)
    epoch = learner.cbstate.history[phase].epochs

    job = addjob!(cb.progress; N = cb.n_steps, 
    description="Epoch: $(epoch)")
    cb.steps_job = job

    update!(cb.epochs_job)
end

function FluxTraining.on(
    event::EpochEnd,
    ::TrainingPhase,
    cb::ProgressCB,
    learner)
    removejob!(cb.progress, cb.steps_job)
end


function FluxTraining.on(
    event::StepEnd,
    ::TrainingPhase,
    cb::ProgressCB,
    learner)
    update!(cb.steps_job)
    render(cb.progress)
end

FluxTraining.stateaccess(::ProgressCB) = (step = (
        loss = Read(),
    ),
    cbstate = (
        history = Read(),
    )    
)
    


# --------------------------- loss plotter callback -------------------------- #

struct LossPlotterCB <: Callback
    training_losses::Vector
    validation_losses::Vector
end

LossPlotterCB() = LossPlotterCB(Float64[], Float64[])

function FluxTraining.on(
    event::EpochEnd,
    phase::TrainingPhase,
    cb::LossPlotterCB,
    learner)
    l =  last(learner.cbstate.metricsepoch[phase][:Loss])[2]
    push!(cb.training_losses, l)
end


function FluxTraining.on(
    event::EpochEnd,
    phase::ValidationPhase,
    cb::LossPlotterCB,
    learner)
    

    l = last(learner.cbstate.metricsepoch[phase][:Loss])[2]
    push!(cb.validation_losses, l)
    length(cb.validation_losses) < 2 && return

    # plot
    y0 = min(minimum(cb.training_losses), minimum(cb.validation_losses))
    y1 = max(maximum(cb.training_losses), maximum(cb.validation_losses))
    T = 1:length(cb.validation_losses)

    l_val = round(last(cb.validation_losses), digits=5)
    l_trn = round(last(cb.training_losses), digits=5)

    plt = lineplot(
        T,
        cb.validation_losses;
        xlabel = "Epoch",
        ylabel = "Loss",
        title = "Loss over epochs",
        name = "Validation ($(l_val))",
        ylim = (y0 - 0.1y0, y1 + 0.1y1),
    ) 
    
    lineplot!(plt, T, cb.training_losses; name = "Training   ($(l_trn))")
    println(plt)
end




FluxTraining.stateaccess(::LossPlotterCB) = (step = (
        loss = Read(),
    ),
    cbstate = (
        metricsepoch = Read(),
    )    
    
)