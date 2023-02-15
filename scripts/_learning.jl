import ForwardDiff
using Term.Progress
import FluxTraining: Callback, EpochEnd, Phase, Read, StepEnd, TrainingPhase, EpochBegin, ValidationPhase
using UnicodePlots



# ---------------------------------------------------------------------------- #
#                                   LEARNING                                   #
# ---------------------------------------------------------------------------- #


function run_model_on_trajectory(model, x_transform, y_transform, warmup, trajectory, can; β=6)
    S, can_decoder = initialize_can_with_warmup(can, warmup, trajectory)    
    x̄, Ω = [], []
    for t in 1:size(trajectory.X, 1)
        x = trajectory.X[t, :]
        ẋ = trajectory.V[t, :]
        J = cover_map_jacobian(x, can)
        v = J*ẋ

        network_input = StatsBase.transform(x_transform, vcat(x, v))
        ω = model(network_input)
        ω = StatsBase.reconstruct(y_transform, ω)

        Ṡ = (can.W * S .+ β * ω .+ b₀) .|> can.σ
        S += (Ṡ - S)/τ

        push!(x̄, can_decoder(S, can)[1])
        push!(Ω, ω)
    end

    return x̄, Ω
end

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
    update!(job; i=-2)
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