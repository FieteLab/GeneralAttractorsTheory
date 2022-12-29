using GeneralAttractors
using GeneralAttractors.Simulations
import IterTools: product

include("../networks/ring.jl")


"""
Run a bunch of simulations varying constant speed v, b₀ and δ
"""
SIMULATE = true

fld_name = "ring_grid_search"
mfld = "ring"
# B = range(4, 6, length = 2) |> collect
B = [1]
D = range(0.05, 0.5, length = 35) |> collect
V = range(0.05, 0.5, length = 4) |> collect
τ = 20


params = product(B, D, V) |> collect
@info "Setting up" length(params)


σ = :softrelu
duration = 150
dt = 0.5
still = 50  # initialization period    
nframes = (Int ∘ round)(duration / dt)


function run_all_sims()
    # ------------------------------- fixed params ------------------------------- #

    count = 0

    for δ in D
        if mfld == "torus"
            can = CAN("torus", cover, n, ξ_t, d_t, k_t; offset_size = δ, σ = σ)
            x₀ = [1, 3.14] # initialize state at center of mfld
            d = map(i -> toruscan.metric(x₀, toruscan.X[:, i]), 1:size(toruscan.X, 2))
            activate = zeros(length(d))
            activate[d.<0.5] .= 1
        
        elseif mfld == "ring"
            can = CAN("ring", cover, n, ξ_r, d_r, k_r; offset_size = δ, σ = σ)
            x₀ = [π/2] # initialize state at position
            d = ringcan.metric.(x₀[1], ringcan.X[1, :])
            activate = zeros(length(d))
            activate[d .< .4] .= 1
        
        else
            error()
        end

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
                    σv = 0,
                    μv = v,
                    vmax = 2v,
                    x₀ = x₀,
                    still = still,
                    modality=:constant,
                )
                simulation = Simulation(can, TJ; η = 0.0, b₀ = b, τ = τ)

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
