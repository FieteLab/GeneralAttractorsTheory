using Plots


using GeneralAttractors
using GeneralAttractors.Simulations



trajectory = Trajectory(torus_attractor; T=250, μ=0.1, θ=0.5)

simulation = Simulation(torus_attractor, trajectory; η = 0.0)

plot(trajectory)


# h = @time run_simulation(
#     simulation;
#     frame_every_n = MODE == :CONSTANT ? 20 : nothing,
#     discard_first_ms = MODE == :CONSTANT ? 0 : 5000,
#     average_over_ms = 20,
# )
