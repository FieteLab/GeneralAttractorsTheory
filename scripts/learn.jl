using Flux
import Flux: mse
using Plots
using FluxTraining
import FluxTraining
using MLDataUtils

using GeneralAttractors
using GeneralAttractors.Simulations
import GeneralAttractors.Simulations: decode_peak_location, Decoder


can_n = 48
warmup_duration = 1000  
trial_duration = 100

# x₀ = nothing # [3.14, 3.14]
b₀ = 1.0
τ = 5.0

α = -270  # scaling factor for ω while generating ground truth data

# ------------------------------ network params ------------------------------ #
d , can_size = 2, can_n^2
n_hidden = 256 
lr = 0.001

n_training_trajectories = 50
n_training_epochs = 250


include("_learning.jl")


# ---------------------------------- get CAN --------------------------------- #
@isdefined(can) || begin
    can = torus_maker(:single; n=can_n)
end



# ------------------------------- generate data ------------------------------ #


"""
    make_data()

Generate a dataset.
It creates a bunch of trajectories and then gets the ground truth 
correct data to train the network on. 
"""
function make_data()
    data = []
    @info "Generating data"
    j = 1
    for i in 1:n_training_trajectories
        # run a "warmup" simulation to initialize the CAN at x₀
        x₀ = rand(can.C.N)
        
        warmup_trajectory = ConstantTrajectory(
            can;
            T=warmup_duration,
            x₀=x₀,
        )

        warmup, _ = run_simulation(
            Simulation(can, warmup_trajectory; η=0.0, b₀=b₀);
            discard_first_ms=warmup_duration-800,
            average_over_ms=0,
        )


        # get data
        trajectory = Trajectory(
            can;
            T=max(trial_duration, 500),  #  needs to be long enough to create a trajectory
            vmax=0.01,
            still=0, 
            x₀=x₀,
        )
        
        _, Ω = generate_groundtruth_data(trajectory, warmup)

        for t in 1:trial_duration
            x_t, v_t = trajectory.X[t, :], trajectory.V[t, :]
            network_input = vcat(x_t, v_t)
            push!(data, (network_input,  Ω[t]))
            j += 1
        end
        println("   $(i)/$(n_training_trajectories)")
    end


    # split train/test 
    train, test = splitobs(shuffleobs(data); at = 0.67)
    @info "Data ready" train test train[1][1] train[2][1]

    return train, test
end

@isdefined(train) || begin 
    train, test = make_data()
end

# --------------------------------- training --------------------------------- #

   
function training_loop(train, test)
    # make network
    model = Chain(
        Dense(2d => n_hidden, relu),  # ẋ, x̂ → ṡ
        Dense(n_hidden => n_hidden, relu),
        Dense(n_hidden => can_size, relu),
    ) 

    # callbacks
    progCB = ProgressCB(length(train), n_training_epochs)

    # train
    trainer = Learner(
                model, mse; optimizer=Flux.ADAM(lr), 
                callbacks=[
                    progCB, 
                    LossPlotterCB(),
                    Recorder(),
                    Metrics(Metric(Flux.mse, phase=TrainingPhase, name="MSE"),) ,
                    throttle(Checkpointer("./data"), EpochEnd, freq=10)
                ], 
                usedefaultcallbacks = false,
                data=(train, test))

    fit!(trainer, n_training_epochs)
    stop!(progCB.progress)
end


encoder = training_loop(train, test)


