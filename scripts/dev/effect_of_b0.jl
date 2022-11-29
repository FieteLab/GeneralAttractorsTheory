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

include("../networks/torus.jl")

"""
Run multiple runs of torus simulations with different b₀
inputs and look at the effect on network bump speed. 
"""

SIMULATE = true
fld_name = "torus_b0_largeker"
dt = 0.5
D = [0.07, 0.1, 0.13] # wegihts offset size
B = range(0.05, 0.5, length=10) # static input
v_actual = 0.025  # actual speed stimulus

# ------------------------------ run simulations ----------------------------- #
if SIMULATE
    dt = 0.5
    duration = 200
    still = 50  # initialization period        

    x₀ = [1, 3.14] # initialize state at center of mfld
    d = map(i -> toruscan.metric(x₀, toruscan.X[:, i]), 1:size(toruscan.X, 2))
    activate = zeros(length(d))
    activate[d.<0.5] .= 1

    # initialize trajectory and simulation
    nframes = (Int ∘ round)(duration / dt)
    for δ in D
        can = CAN(
            "torus",
            cover,
            n,
            ξ_t,
            d_t,
            k_t;
            offset_size = δ,
            # Ω = Ω
        )
        trajectory = Trajectory(
            can;
            T = nframes,
            dt = dt,
            σv = v_actual,
            μv = v_actual,
            vmax = v_actual,
            σθ = 0.0,
            θ₀ = 0,
            x₀ = 1, y₀ = 1,
            still = still,
        )

        for (i, b₀) in enumerate(B)
            println(Panel("Running sim $i/$(length(B))", style = "red", justify = :center))

            simulation = Simulation(can, trajectory; η = 0.0, b₀ = b₀)

            # run
            h, X̄ = @time run_simulation(
                simulation;
                frame_every_n = nothing,
                discard_first_ms = 50,
                average_over_ms = 1,
                fps = 10,
                s₀ = 1.0 .* activate,
                savefolder = fld_name,
                savename = "δ_$(δ)_run_$(i)",
            );
            # println(size(h.S))
            plot_trajectory_and_decoded(trajectory, X̄) |> display
        end
    end
end


# ------------------------------- run analysis ------------------------------- #

plt = plot(xlabel="b0", ylabel="mfld displacement", title="mfld displacement")
plt2 = plot(xlabel="b0", ylabel="on mfld speed", )

function ∑(x)
    n, _, m = size(x)
    real.(reshape(mean(x, dims = 2), (n, m)))
end

for (δ, color) in zip(D, (:black, :red, :green))
    V, S = [], []  # will store average state speed for each simulation

    for (i, b₀) in enumerate(B)
        X = load_data(fld_name, "δ_$(δ)_run_$(i)_torus_sim_decoded_X")

        dx = mean(diff(X[:, 1]))
        push!(V, dx/δ)

        # load state history
        history = load_simulation_history(fld_name, "δ_$(δ)_run_$(i)_torus_history")
        s = ∑(history.S)[:, 10:end]

        # get peak location speed
        peak_location = hcat(
            map(
                st -> decode_peak_location(st, toruscan), 
                eachcol(s)
                )...
            )

        on_mfld_speed = map(
            i -> toruscan.metric(
                peak_location[:, i], peak_location[:, i-1]
                ), 
            2:size(peak_location, 2)
        )
        average_speed = sum(on_mfld_speed)/length(on_mfld_speed)
        plot(on_mfld_speed, ylim=[0, 5]) |> display
        @info "cacca" sum(on_mfld_speed)

        push!(S, average_speed)
        # push!(S, mean(s))

    end

    plot!(plt, B, V, lw=1, color=color, label="offset: $δ")
    scatter!(plt, B, V, ms=3, color="white", msc=color, label=nothing)

    plot!(plt2, B, S, lw=1, color=color, label="offset: $δ")
    scatter!(plt2, B, S, ms=3, color="white", msc=color, label=nothing)
end


hline!(plt2, [v_actual], lw=3, alpha=.5, color=:blue, label="V actual")
plot(plt2,  size=(800, 600))