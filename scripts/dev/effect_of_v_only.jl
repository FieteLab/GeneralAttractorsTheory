using Plots


using GeneralAttractors
using GeneralAttractors.Simulations

using Statistics
using Term
using Distances
using GeneralAttractors.ManifoldUtils
using GeneralAttractors.Analysis
import GeneralAttractors.Analysis.ManifoldAnalysis: population_average
import GeneralAttractors.Simulations: decode_peak_location, plot_trajectory_and_decoded
import GeneralAttractors.Analysis: get_bump_speed

include("../networks/torus.jl")

"""
Run multiple runs of torus simulations with different b₀
and input velocities to see how it affects the bump speed
"""

SIMULATE = true
fld_name = "torus_v_only"
can = toruscan

dt = 0.5
V = range(0.001, 0.05, length = 7) # speed stimuli


# ------------------------------ run simulations ----------------------------- #
if SIMULATE
    dt = 0.5
    duration = 250
    still = 50  # initialization period        

    x₀ = [1, 3.14] # initialize state at center of mfld
    d = map(i -> toruscan.metric(x₀, toruscan.X[:, i]), 1:size(toruscan.X, 2))
    activate = zeros(length(d))
    activate[d.<0.5] .= 1

    # initialize trajectory and simulation
    nframes = (Int ∘ round)(duration / dt)
    for (j, v) in enumerate(V)
        "\n"^6 / hLine() |> print
        println(Panel("Running sim $j/$(length(V))", style = "red", justify = :center))

        TJ = Trajectory(
            can;
            T = nframes,
            dt = dt,
            μv = v,
            vmax = v,
            still = still,
            modality=:constant
        )


        simulation = Simulation(can, TJ; η = 0.0, b₀ = 7.0, τ=5)

        # run
        h, X̄ = @time run_simulation(
            simulation;
            frame_every_n = nothing,
            discard_first_ms = still,
            average_over_ms = 1,
            fps = 10,
            s₀ = 1.0 .* activate,
            savefolder = fld_name,
            savename = "v_$(v)",
        )

    end
end



S = []
trajectory_V = []
for v in V
    sim_name = "v_$(v)_torus"
    push!(trajectory_V, v)

    # load state history
    push!(S, get_bump_speed(can, fld_name, sim_name*"_history"))
end


plt = plot(
    xlabel = "velocity input",
    ylabel = "bump  velocity",
    title = "δ constant, b₀ constant",
)

plot!(plt, V, S, lw = 2, color = :black, label = nothing)
scatter!(plt, V, S, ms = 5, color = "white", msc = :black, label = nothing)

plot!(plt, V, V, lw = 2, color = :black, alpha = 0.5, ls = :dash, label = nothing)
plot(plt, size = (800, 600), grid = false, aspect_ratio = :equal)
