using Flux
import Flux: mse
using Plots
using FluxTraining
import FluxTraining
using MLDataUtils
using JLD2

using GeneralAttractors
using GeneralAttractors.Simulations
import GeneralAttractors.Simulations: decode_peak_location, Decoder


include("_learning.jl")


can_n = 48
warmup_duration = 1000  
trial_duration = 1500
x₀= [3.14, 3.14]
    
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


trajectory = Trajectory(
    can;
    T=max(trial_duration, 500),  #  needs to be long enough to create a trajectory
    vmax=0.01,
    still=0, 
    x₀=x₀,
)

# get ground truth data
println("Generating ground truth data")
x̄, Ω = generate_groundtruth_data(trajectory, warmup)

# get data through network
model = loadmodel("data/checkpoint_epoch_056_loss_0.0039244010565047905.bson")
x_transform = load_object("./data/x_transform.jld2")
y_transform = load_object("./data/y_transform.jld2")

S, can_decoder = initialize_can(warmup, trajectory)    
println("Running model on trajectory")
model_x̄, model_Ω = run_model_on_trajectory(
    model, x_transform,  y_transform, warmup, trajectory, can
)


plot(eachcol(trajectory.X)..., lw=4, color=:black, label="trajectory")
plot!(eachrow(hcat(x̄...))..., lw=2, color=:red, label="ground truth")
plot!(eachrow(hcat(model_x̄...))..., lw=2, color=:blue, label="model")