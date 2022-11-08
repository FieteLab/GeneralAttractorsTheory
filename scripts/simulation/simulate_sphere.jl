using GeneralAttractors
using Plots
using GeneralAttractors.Simulations


dt = 0.5
duration = 800

nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(sphere_attractor; T = nframes)
simulation = Simulation(sphere_attractor, trajectory)

h = run_simulation(
    simulation,
    frame_every_n = 20,
    discard_first_ms = 0,
    average_over_ms = 20,
)
