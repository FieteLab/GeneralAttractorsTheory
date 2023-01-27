using Plots


using GeneralAttractors
using GeneralAttractors.Simulations
import GeneralAttractors.Simulations: plot_trajectory_and_decoded
import GeneralAttractors: animate_simulation_data

include("settings.jl")


"""
TODO: 
- [ ] make all CANs work on the same traj
- [ ] save a simulation data for all cans
- [ ] save an animation of the trajectory, the decoded trajectory, the bump and the neuron state space
"""


# pyplot()

# ---------------------------------- get CAN --------------------------------- #
dt = 0.5
duration = 1000
still = 50  # initialization period        

network = "torus"

can, x₀_net, x₀_traj, embedding = if network == "torus"
    can = torus_maker(:defult; n=48, α=31)
    x₀_net, x₀_traj = [3.14, 3.14], [-10, 0]
    can, x₀_net, x₀_traj, torus_embedding
elseif network == "cylinder"
    can = cylinder_maker(:default; n=48, α=30)
    x₀_net, x₀_traj = [3.14, 0.0] , [-10, 0]
    can, x₀_net, x₀_traj, cylinder_embedding
elseif network == "plane"
    can = plane_maker(:default; n=48, α=140, offset_size=0.05)
    x₀_net, x₀_traj = [0.0, 0.0], [-10, 0]
    can, x₀_net, x₀_traj, plane_embedding
end

d = map(i -> can.metric(x₀_net, can.X[:, i]), 1:size(can.X, 2))
activate = zeros(length(d))
activate[d.<1.5] .= 1


# ------------------------ make simulation trajecotry ------------------------ #
nframes = (Int ∘ round)(duration / dt)
v_mag = (cos.(range(0, 2π - .1, length=nframes)) ./ 2 .+ .5)
vx = sin.(range(0, 1, length=nframes)) .* v_mag
# vx = collect(1:nframes) ./ nframes .* v_mag
vy = cos.(range(0, 2π - .1, length=nframes)) .* v_mag

trajectory = Trajectory(
    can;
    T = nframes,
    dt = dt,
    σv = 5,
    # μv = [0.1, 0.07],
    μv = 0,
    vmax = 0.05,
    still = still,
    x₀ = x₀_traj,
    smoothing_window = 501,
    Vs = [vx, vy],
)
plot(trajectory) |> display
simulation = Simulation(can, trajectory; η = 0.0, b₀ = 1.0);

# --------------------------------- simulate --------------------------------- #
h, X̄ = @time run_simulation(
    simulation;
    discard_first_ms = still,
    average_over_ms = 0,
    s₀ = 1.0 .* activate,
    # savefolder = "abstract",
    # savename = "torus_$i",
)


# --------------------------------- visualie --------------------------------- #
plot_trajectory_and_decoded(trajectory, X̄) |> display
animate_simulation_data(can, trajectory, h, X̄, embedding, 
        (supervisor.projectdir / "plots" /"$(network)_sim_traj.gif").path
)