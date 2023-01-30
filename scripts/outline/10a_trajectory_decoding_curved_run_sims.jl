using Plots


using GeneralAttractors
using GeneralAttractors.Simulations
using GeneralAttractors.Kernels
import GeneralAttractors.Simulations: plot_trajectory_and_decoded
import GeneralAttractors: animate_simulation_data

include("settings.jl")
move_to_datadir(supervisor, "path_int")
tag = "decoding_data"

# ---------------------------------- get CAN --------------------------------- #

duration = 500
still = 50  # initialization period        

network = "mobius"
can, x₀_traj, embedding = make_path_int_can(network)


# ------------------------ make simulation trajecotry ------------------------ #
nframes = (Int ∘ round)(duration / dt)
v_mag = (cos.(range(0, 2π - .1, length=nframes)) ./ 2 .+ .5)
vx = sin.(range(0, 1, length=nframes)) .* v_mag 
vy = cos.(range(0, 2π - .1, length=nframes)) .* v_mag 

trajectory = Trajectory(
    can;
    T = nframes,
    dt = dt,
    σv = 5,
    μv = 0,
    vmax = 0.05,
    still = still,
    x₀ = x₀_traj,
    smoothing_window = 501,
    Vs = [vx, vy],
)
trajplot = plot(trajectory)
save_plot(supervisor, trajplot, "10_path_int_$(network)_traj")
simulation = Simulation(can, trajectory; η = 0.0, b₀ = 1.0);


# --------------------------------- simulate --------------------------------- #

# activate = get_can_initialization_weights(trajectory, can)

# h, X̄ = @time run_simulation(
#     simulation;
#     discard_first_ms = still,
#     average_over_ms = 0,
#     s₀ = 1.0 .* activate,
# )


# meta = metadata = Dict(
#     "network" => network,
#     "dt" => dt,
#     "duration" => duration,
#     "still" => still,
#     "tag" => tag,
# )


# store_data(supervisor, "simulations"; fmt = "jld2", name="$(network)_sim_data", data = Dict("h" => h, "trajectory"=>trajectory), metadata=meta)

# # --------------------------------- visualie --------------------------------- #
# plot_trajectory_and_decoded(trajectory, X̄) |> display
# # animate_simulation_data(can, trajectory, h, X̄, embedding, 
# #         (supervisor.projectdir / "plots" /"path_int_$(network)_sim.gif").path
# # )