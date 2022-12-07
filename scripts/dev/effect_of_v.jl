using Plots


using GeneralAttractors
using GeneralAttractors.Simulations

using LinearAlgebra: norm
using Statistics
using Term
using Distances
using GeneralAttractors.ManifoldUtils
using GeneralAttractors.Analysis
import GeneralAttractors.Analysis.ManifoldAnalysis: population_average, tda_on_pointcloud
import GeneralAttractors.Simulations: decode_peak_location, plot_trajectory_and_decoded

import MyterialColors: indigo, salmon_dark, Palette

include("../networks/torus.jl")

"""
Run multiple runs of torus simulations with different b₀
and input velocities to see how it affects the bump speed
"""

SIMULATE = false
fld_name = "torus_v"
dt = 0.5
V = range(0.01, 0.1, length=4) # speed stimuli
B = range(0.05, 0.5, length=8) # static input

colors = getfield.(Palette(indigo, salmon_dark; N=length(V)).colors, :string)

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
        can = CAN(
            "torus",
            cover,
            n,
            ξ_t,
            d_t,
            k_t;
            offset_size = 0.2,
            # Ω = Ω
        )
        TJ = Trajectory(
            can;
            T = nframes,
            dt = dt,
            σv = v,
            μv = v,
            vmax = v,
            σθ = 0.0,
            θ₀ = 0,
            x₀ = 1, y₀ = 1,
            still = still,
        )

        for (i, b₀) in enumerate(B)
            println(Panel("Running sim $i/$(length(B)) | CAN: $j/$(length(V))", style = "red", justify = :center))

            simulation = Simulation(can, TJ; η = 0.0, b₀ = b₀)

            # run
            h, X̄ = @time run_simulation(
                simulation;
                frame_every_n = nothing,
                discard_first_ms = 50,
                average_over_ms = 1,
                fps = 10,
                s₀ = 1.0 .* activate,
                savefolder = fld_name,
                savename = "v_$(v)_run_$(i)",
            );

            # plot_trajectory_and_decoded(TJ, X̄) |> display
        end
    end
end


# ------------------------------- run analysis ------------------------------- #


can = CAN(
    "torus",
    cover,
    n,
    ξ_t,
    d_t,
    k_t;
    offset_size = 0.1,
    )


plt = plot(xlabel="b0", ylabel="on mfld speed", title="δ constant, v varies")

function ∑(x)
    n, _, m = size(x)
    real.(reshape(mean(x, dims = 2), (n, m)))
end

for (v, color) in zip(V, colors)
    S = []

    for (i, b₀) in enumerate(B)
        # load state history
        history = load_simulation_history(fld_name, "v_$(v)_run_$(i)_torus_history")
        s = ∑(history.S)[:, 10:end]

        # get peak location speed
        peak_location = hcat(
            map(
                st -> decode_peak_location(st, can), 
                eachcol(s)
                )...
            )

        on_mfld_speed = map(
            i -> can.metric(
                peak_location[:, i], peak_location[:, i-1]
                ), 
            2:size(peak_location, 2)
        )
        average_speed = sum(on_mfld_speed)/(length(on_mfld_speed)*dt)  # tot displacement over time

        push!(S, average_speed)
    end

    plot!(plt, B, S, lw=2, color=color, label="actual v: $v")
    scatter!(plt, B, S, ms=5, color="white", msc=color, label=nothing)
    hline!(plt, [v], lw=1, color=color, label=nothing, linestyle=:dash, alpha=.5)
end


plot(plt,  size=(800, 600), grid=false)