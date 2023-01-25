import GeneralAttractors.Can: SingleCAN
using Ripserer
import OrderedCollections: OrderedDict
using Term.Tables
import GeneralAttractors.Analysis.ManifoldAnalysis: 
        tda_on_pointcloud, estimate_intrinsic_dimensionality


function make_standard_single_torus_can(; kernel_name = :DoE)
    can_maker = network_makers["torus"]
    kernel_params = Dict(
        k => mean(v) for (k, v) in kernels_parameters_range["torus"][kernel_name]
    )
    kernel = kernels[kernel_name](; kernel_params...)
    return can_maker(:single; k=kernel)
end


"""
Get a random initial condition for a CAN
"""
function random_init(can; x₀=nothing)
    x₀ = something(x₀, rand(can.C.N))
    d = map(i -> can.metric(x₀, can.X[:, i]), 1:size(can.X, 2))
    activate = zeros(length(d))
    activate[d.<0.5] .= 1
    return x₀, activate
end

"""
    constant_traj_si(can, duration, dt, still, τ, b₀)

Get a sim for a single CAN with a constant trajectory.
"""
function constant_traj_sim(can::SingleCAN, duration, dt, still, τ, b₀; η=0.1)
    # initialize trajectory and simulation
    nframes = (Int ∘ round)(duration / dt)
    trajectory = ConstantTrajectory(
        can;
        T = nframes,
        still = still,
    )

    return Simulation(can, trajectory; b₀ = b₀, τ = τ, η=0.1)
end

"""
    simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀)

Get and run a simulation for a SingleCAN wwith random initial condition
and constant input trajectory.
"""
function simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀; x₀=nothing, 
    η=0.0, kwargs...
    )
    x₀, activate = random_init(can; x₀=x₀)
    sim = constant_traj_sim(can, duration, dt, still, τ, b₀; η=η)
    h, X = run_simulation(    
        sim;
        discard_first_ms = still,
        average_over_ms = 1,
        s₀ = 1.0 .* activate,
        kwargs...
    );
    return h, X
end



# --------------------------------- topology --------------------------------- #

"""
Load a subset of the data from the supervisor with multiple simulations and concatenate
activations. 
"""
function load_and_concat_activations(; filters...)
    _, data = ProjectSupervisor.fetch(supervisor; filters...)
    
    # stack activations over time
    X = hcat(
        map(
            d -> d["S"][:, 1, end-5:end], data
        )...
    )
    @info "Loaded $(length(data)) simulations." X
    return X
end


"""
do PCA dimensionality reduction
"""
do_pca(X, params) = Analysis.pca_dimensionality_reduction(X, params)[2]

"""
do Isomap dimensionality reduction
"""
do_isomap(X, params) = isomap_dimensionality_reduction(X, params)[2]


# ------------------------------------ tda ----------------------------------- #
"""
Get the number of persistent (not noise) features 
from a TDA persistence diagram based on the 
largest itnerval approach.
"""
function get_n_persistence_features(intervals)
    per = persistence.(intervals) |> sort
    gaps = per[2:end] .- per[1:end-1]
    largest = argmax(gaps)
    return length(per) - largest
end

"""
    do_tda(
        supervisor::Supervisor, 
        data_filters::Dict,
        save_plot_name::String;
        max_d = 1
        )

Load pointcloud data (generally after isomap dim red)
and run TDA on it. save barcode and persistence diagram plot
and give the number of relevant features. 
"""
function do_tda(
    supervisor::Supervisor, 
    data_filters::Dict,
    save_plot_name::String;
    max_d = 1,
    tresh = 20
    )
    tda_params = AnalysisParameters(
        tda_threshold = tresh,       # threshold to reduce TDA computation 
        tda_downsample_factor = 10,        # temporal downsampling of data for TDA
        tda_dim_max = max_d,        # max feature dimension, starting at 0
    )
    
    # load data
    _, M = ProjectSupervisor.fetch(supervisor; data_filters...) 
    @assert length(M) == 1 length(M)
    M = M[1]
    # @assert size(M, 1) == 10 size(M)

    # run TDA
    tda_colors = [green_dark, indigo, salmon]
    tda_model, tda_barcode_plot = tda_on_pointcloud(M, tda_params; 
            ms=5, color = tda_colors[1:max_d+1],
            plot_font_size_kwargs...)
    save_plot(supervisor, tda_barcode_plot,save_plot_name);


    # get the number of persistence features based on the largest gap in the persitence
    features = OrderedDict(
        map(d -> d=>get_n_persistence_features(tda_model[d]), 1:max_d+1)...
    )

    # print the result
    Panel(
        Table(Dict(
        "dim" => collect(keys(features)),
        "n_features" => collect(values(features)),
    ));
        title = data_filters[:can],
        title_style = "blue bold",
        style="bright_yellow",
        fit=false,
        width = 50
    ) |> print
end