using DrWatson
@quickactivate "GeneralAttractors"


using GeneralAttractors
using GeneralAttractors.Simulations
using Term
install_term_stacktrace(hide_frames = true)

using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp, moving_average
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: plot_trajectory_and_decoded, plot_on_mfld_trajectory_and_history
import GeneralAttractors.Analysis: get_bump_speed

"""
Run a large number of simulations for a single copy CAN with different initial conditions.
Save the data to make sure  
"""

# include("../networks/torus.jl")


# --------------------------------- functions -------------------------------- #
function random_init(can)
    x₀ = can.C.N |> rand
    d = map(i -> can.metric(x₀, can.X[:, i]), 1:size(can.X, 2))
    activate = zeros(length(d))
    activate[d.<0.5] .= 1
    return x₀, activate
end

function setup(can, duration, dt, still, τ, b₀, x₀)
    # initialize trajectory and simulation
    nframes = (Int ∘ round)(duration / dt)
    trajectory = Trajectory(
        can;
        T = nframes,
        dt = dt,
        x₀ = x₀,
        μv = 0,
        still = still,
        modality = :constant
    )

    return Simulation(can, trajectory; η = 0.0, b₀ = b₀, τ = τ)
end

function run(simulation, still, activate)
    h, X̄ = @time run_simulation(    
        simulation;
        frame_every_n = nothing,
        discard_first_ms = still,
        average_over_ms = 1, 
        fps = 4,
        s₀ = 1.0 .* activate,
        savefolder = "torus",
        savename = "test",
    );
end

# ---------------------------------------------------------------------------- #
#                                      RUN                                     #
# ---------------------------------------------------------------------------- #
# select model


params = Dict{Symbol, Any}(
    :can => toruscan,
    :N_sims => 10,
    :dt => 0.5,
    :duration => 200,
    :still => 100,
    :τ => 5.0,
    :b₀ => 0.5  ,
)

for i in 1:params[:N_sims]
    x₀, activate = random_init(params[:can])
    # sim = setup(params[:can], params[:duration], params[:dt], params[:still], params[:τ], params[:b₀], x₀)
    # h, X̄ = run(sim, params[:still], activate)

    save_info = copy(params)
    save_info[:x₀] = x₀
    save_info[:i] = i


    save_path = datadir("simulations", "random_inits", savename(save_info, "npz"))
    wsave(
        save_path,
        zeros(100),   
    )
end