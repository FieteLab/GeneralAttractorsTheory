using Plots


using GeneralAttractors
using GeneralAttractors.Simulations

using Statistics
using Term
using Distances
using GeneralAttractors.ManifoldUtils
using GeneralAttractors.Analysis
import GeneralAttractors.Analysis.ManifoldAnalysis: population_average, tda_on_pointcloud

include("../networks/torus.jl")

"""
Run multiple runs of torus simulations with different b₀
inputs and look at the effect on network bump speed. 
"""

SIMULATE = false
D = [0.075, 0.1, 0.125] # wegihts offset size
B = range(0.01, 1, length=20) # static input
v_actual = 0.2  # actual speed stimulus

# ------------------------------ run simulations ----------------------------- #
if SIMULATE
    dt = 0.5
    duration = 700
    still = 50  # initialization period        


    x₀ = [3.14, 3.14] # initialize state at center of mfld
    d = map(i -> toruscan.metric(x₀, toruscan.X[:, i]), 1:size(toruscan.X, 2))
    activate = zeros(length(d))
    activate[d.<0.5] .= 1

    # initialize trajectory and simulation
    nframes = (Int ∘ round)(duration / dt)
    trajectory = Trajectory(
        toruscan;
        T = nframes,
        dt = dt,
        σv = v_actual,
        μv = v_actual,
        vmax = v_actual,
        σθ = 0.0,
        θ₀ = 0,
        still = still,
    )

    for δ in D
        toruscan = CAN(
            "torus",
            cover,
            n,
            ξ_t,
            d_t,
            k_t;
            offset_size = δ,
            # Ω = Ω
        )
        for (i, b₀) in enumerate(B)
            println(Panel("Running sim $i/$(length(B))", style = "red", justify = :center))

            simulation = Simulation(toruscan, trajectory; η = 0.0, b₀ = b₀)

            # run
            h, X̄ = @time run_simulation(
                simulation;
                frame_every_n = nothing,
                discard_first_ms = 100,
                average_over_ms = 10,
                fps = 10,
                s₀ = 1.0 .* activate,
                savefolder = "torus_b0",
                savename = "δ_$(δ)_run_$(i)",
            );
        end
    end
end


# ------------------------------- run analysis ------------------------------- #

plt = plot(xlabel="b0", ylabel="mfld displacement", title="mfld displacement", )
plt2 = plot(xlabel="b0", ylabel="mean activity", title="mean activity", )

for (δ, color) in zip(D, (:black, :red, :green))
    V,S = [], []  # will store average state speed for each simulation

    for (i, b₀) in enumerate(B)
        X = load_data("torus_b0", "δ_$(δ)_run_$(i)_torus_sim_decoded_X")
        dx = sum(diff(X[:, 1]))
        push!(V, dx/δ)

        
        history = load_simulation_history("torus_b0", "δ_$(δ)_run_$(i)_torus_history")
        s = real.(population_average(history; skip = 20))

        push!(S, mean(s))
    end

    plot!(plt, B, V, lw=1, color=color, label="offset: $δ")
    scatter!(plt, B, V, ms=3, color="white", msc=color, label=nothing)

    plot!(plt2, B, S, lw=1, color=color, label="offset: $δ")
    scatter!(plt2, B, S, ms=3, color="white", msc=color, label=nothing)
end

plot(plt, plt2, layout=(2, 1), size=(800, 600))