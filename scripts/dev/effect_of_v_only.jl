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


include("../networks/torus.jl")

"""
Run multiple runs of torus simulations with different b₀
and input velocities to see how it affects the bump speed
"""

SIMULATE = true
fld_name = "torus_v_only"
dt = 0.5
V = range(0.001, 1.0, length=5) # speed stimuli


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


        can = CAN(
            "torus",
            cover,
            n,
            ξ_t,
            d_t,
            k_t;
            offset_size = 1.5,
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


        simulation = Simulation(can, TJ; η = 0.0, b₀ = 7.0)

        # run
        h, X̄ = @time run_simulation(
            simulation;
            frame_every_n = nothing,
            discard_first_ms = 50,
            average_over_ms = 1,
            fps = 10,
            s₀ = 1.0 .* activate,
            savefolder = fld_name,
            savename = "v_$(v)",
        );

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



function ∑(x)
    n, _, m = size(x)
    real.(reshape(mean(x, dims = 2), (n, m)))
end


S = []
trajectory_V = []
for v in V
    # get trjectory's speed
    X = load_data(fld_name, "v_$(v)_torus_sim_trajectory_X")
    traj_v = sum(diff(X[:, 1]))/(size(X, 1)*dt)
    push!(trajectory_V, traj_v)

    # load state history
    history = load_simulation_history(fld_name, "v_$(v)_torus_history")
    s = ∑(history.S)[:, 30:end]

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
    push!(S, average_speed/(0.5))
end


plt = plot(xlabel="velocity input", ylabel="bump  velocity", title="δ constant, b₀ constant")

plot!(plt, V, S, lw=2, color=:black, label=nothing)
scatter!(plt, V, S, ms=5, color="white", msc=:black, label=nothing)

plot!(plt, V, V, lw=2, color=:black, alpha=.5, ls=:dash, label=nothing)
plot(plt,  size=(800, 600), grid=false, aspect_ratio=:equal)