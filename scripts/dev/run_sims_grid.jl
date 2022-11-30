using GeneralAttractors
using GeneralAttractors.Simulations
import IterTools: product

include("../networks/torus.jl")


"""
Run a bunch of simulations varying constant speed (v), b₀ and δ
"""
SIMULATE = true

B = range(0.1, 0.4, length=4)
D = range(0.1, 0.7, length=7)
V = range(0.1, 1.0, length=19)
params = product(B, D, V) |> collect
@info "Setting up" length(params)


function run_all_sims()
    fld_name = "params_sims_lin"

    # ------------------------------- fixed params ------------------------------- #
    dt = 0.5
    duration = 250
    still = 50  # initialization period    
    nframes = (Int ∘ round)(duration / dt) 

    x₀ = [1, 3.14] # initialize state at center of mfld
    d = map(i -> toruscan.metric(x₀, toruscan.X[:, i]), 1:size(toruscan.X, 2))
    activate = zeros(length(d))
    activate[d.<0.5] .= 1

    count = 0

    for δ in D
        can = CAN(
            "torus",
            cover,
            n,
            ξ_t,
            d_t,
            k_t;
            offset_size = δ,
        )

        for b in B
            for v in V
                count += 1
                println(Panel("Running sim $count/$(length(params))", style = "red", justify = :center))

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
                simulation = Simulation(can, TJ; η = 0.0, b₀ = b)

                # run
                h, X̄ = @time run_simulation(
                    simulation;
                    frame_every_n = nothing,
                    discard_first_ms = 50,
                    average_over_ms = 1,
                    fps = 10,
                    s₀ = 1.0 .* activate,
                    savefolder = fld_name,
                    savename = "v_$(v)_δ_$(δ)_b_$b",
                );
            end
        end
    end
end


SIMULATE && run_all_sims()