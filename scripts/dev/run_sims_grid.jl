using GeneralAttractors
using GeneralAttractors.Simulations
import IterTools: product

include("../networks/torus.jl")


"""
Run a bunch of simulations varying constant speed v, b₀ and δ
"""
SIMULATE = true

fld_name = "params_grid_search_tanh"
B = range(4, 6, length = 2) |> collect
D = range(0.8, 1.4, length = 6) |> collect
V = range(0.1, 0.7, length = 20) |> collect


params = product(B, D, V) |> collect
@info "Setting up" length(params)

σ = :tanh
duration = 150
dt = 0.5
still = 50  # initialization period    
nframes = (Int ∘ round)(duration / dt)


function run_all_sims()
    # ------------------------------- fixed params ------------------------------- #
    x₀ = [1, 3.14] # initialize state at center of mfld
    d = map(i -> toruscan.metric(x₀, toruscan.X[:, i]), 1:size(toruscan.X, 2))
    activate = zeros(length(d))
    activate[d.<0.5] .= 1

    count = 0

    for δ in D
        can = CAN("torus", cover, n, ξ_t, d_t, k_t; offset_size = δ, σ = σ)

        for b in B
            for v in V
                count += 1
                println(
                    Panel(
                        "Running sim $count/$(length(params))",
                        style = "red",
                        justify = :center,
                    ),
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
                    x₀ = [1, 1],
                    still = still,
                )
                simulation = Simulation(can, TJ; η = 0.0, b₀ = b, τ = 5)

                # run
                h, X̄ = @time run_simulation(
                    simulation;
                    frame_every_n = nothing,
                    discard_first_ms = still,
                    average_over_ms = 0,
                    fps = 10,
                    s₀ = 1.0 .* activate,
                    savefolder = fld_name,
                    savename = "v_$(v)_δ_$(δ)_b_$b",
                )
            end
        end
    end
end


SIMULATE && run_all_sims()
