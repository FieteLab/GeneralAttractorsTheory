import GeneralAttractors.Can: SingleCAN
using Ripserer
import OrderedCollections: OrderedDict
using Term.Tables
import GeneralAttractors.Analysis.ManifoldAnalysis: 
        tda_on_pointcloud, estimate_intrinsic_dimensionality
using MultivariateStats, ManifoldLearning

# ------------------------------ networks makers ----------------------------- #

function make_single_can(network::String; kwargs...)
    n = network ∈ ("ring", "line") ? 256 : 48
    return network_makers[network](:single; n=n, kwargs...)
end

"""
Generate a CAN with multiple copies for path integration stuff.
"""
function make_path_int_can(network; funky=false, random_x0=false, kwargs...) 
    return if network == "line"
        can = line_maker(:default; n=256, α=21, kwargs...)
        x₀_traj = random_x0 ? nothing : [2]
        can, x₀_traj, line_embedding
    elseif network == "ring"
        can = ring_maker(:default; n=256, α=21, kwargs...)
        x₀_traj = random_x0 ? nothing : [π/2]
        can, x₀_traj, ring_embedding
    
    elseif network == "torus"
        can = if funky
                torus_maker(:defult; n=48, α=75, offset_size=0.15, use_offset_fields=true, kwargs...) # for curly vfield α=325
            else
                torus_maker(:defult; n=48, α=18, offset_size=0.3, use_offset_fields=false, kwargs...)
        end
        x₀_traj = random_x0 ? nothing : [-20, -15]
        can, x₀_traj, torus_embedding

    elseif network == "cylinder"
        can = cylinder_maker(:default; n=48, α=30, kwargs...)
        x₀_traj = random_x0 ? nothing : [-20, -15]
        can, x₀_traj, cylinder_embedding

    elseif network == "plane"
        can = plane_maker(:default; 
                n=48, α=65, offset_size=0.1,
                )
        x₀_traj = random_x0 ? nothing : [-20, -15]
        can, x₀_traj, plane_embedding

    elseif network == "mobius"
        can = mobius_maker(:defult; n=48, α=130, offset_size=0.2, kwargs...)
        x₀_traj = random_x0 ? nothing : [0.5, 0.5]
        can, x₀_traj, mobius_embedding

    elseif network == "sphere"
        can = sphere_maker(:default; n=48, α=235, use_offset_fields=funky, kwargs...) # 225
        x₀_traj = random_x0 ? nothing : ([0, 1, 0] ./ norm([0, 1, 0]))
        can, x₀_traj, sphere_embedding
    elseif network == "klein"
        can = klein_maker(:default; n=48, α=235, use_offset_fields=funky, kwargs...) # 225
        x₀_traj = random_x0 ? nothing : [0.5, 0.5]
        can, x₀_traj, klein_embedding
    end
end




# ----------------------------------- misc ----------------------------------- #

"""
Create mask to initialize CAN at.
"""
function get_can_initialization_weights(trajectory, can; δ=1.5)
    x₀_net = trajectory.X̄[1, :]
    d = map(i -> can.metric(x₀_net, can.X[:, i]), 1:size(can.X, 2))
    activate = zeros(length(d))
    activate[d .< δ] .= 1
    return activate
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

    return Simulation(can, trajectory; b₀ = b₀, τ = τ, η=η)
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
function load_and_concat_activations(; expected_n=nothing, filters...)
    _, data = ProjectSupervisor.fetch(supervisor; filters...)
    println(length(data))

    # if expected_n !== nothing
    #     @assert length(data) == expected_n "Expected $(expected_n) simulations, got $(length(data))."
    # end
    
    # stack activations over time
    data = filter(d -> "S" in keys(d), data)
    X = hcat(
        map(
            d -> d["S"][:, 1, end-5:end], data
        )...
    )
    @info "Loaded $(length(data)) simulations."
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
    # gaps = per[2:end] .- per[1:end-1]
    # largest = argmax(gaps)
    # return length(per) - largest

    # get the number of internvals > 1
    return length(per[per .> 6])
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
    network::String,
    η::Float64,
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

    
    @assert length(M) >= 1 "Loaded too much data ($(length(M)) entries) with filters $data_filters"
    @assert length(M) > 0 "Loaded no data with filters $data_filters"
    M = M[1] isa Matrix ? M[1] : M[2]
    @assert M isa Matrix "Expected a matrix of data, got $(typeof(M))."

    # @assert size(M, 1) == 10 size(M)

    # run TDA
    tda_colors = [green_dark, indigo, salmon]
    tda_model, tda_barcode_plot = tda_on_pointcloud(M, tda_params; 
            ms=5, color = tda_colors[1:max_d+1],
            plot_font_size_kwargs...)



    # get the number of persistence features based on the largest gap in the persitence
    features = nothing
    try
        features = OrderedDict(
            map(d -> d=>get_n_persistence_features(tda_model[d]), 1:max_d+1)...
        )
    catch
        @warn "TDA failed for $(data_filters[:can])."
        return
    end

    # save persistence plot
    plt = plot(tda_barcode_plot, size=(800, 600), plot_title="$(network) - Betti: " * string(values(features)))
    save_plot(supervisor, plt, save_plot_name);
    display(plt)

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


# ---------------------------------------------------------------------------- #
#                                   PLOTTING                                   #
# ---------------------------------------------------------------------------- #


"""
Make a frame of a gif to show path integration dynamics:
i. decoded trajectory
ii. activation on neural lattice
iii. activation in ISO embedding of state space
"""
function path_int_gif_frame(d, trajectory, X, coord3d, S, can, w_x, w_y, M, S_embedd, fnum, skipframes, hist_frame; p3_extent=0.5, spatial_downsampling=1)
    # plot trajectory & decoded
    p1 = plot(trajectory, fnum+skipframes)
    plot!(
        eachrow(X[:, 1:fnum])...,
        lw = 6, color=:red, alpha=.5,
        label = "decoded",
        xlabel = "m₁", ylabel = "m₂",
        )
    s =  S[:, hist_frame]

    # plot activation on neural lattice
    p2 = if d == 2
        contourf(
            w_x, w_y, Matrix(reshape(s, can.n)'), 
            xlabel = "θ₁", ylabel = "θ₂",
            aspect_ratio = :equal,
            linewidth = 0.25,
            xlim = (minimum(w_x)-0.5, maximum(w_x)+0.5), ylim = (minimum(w_y)-0.5, maximum(w_y)+0.5),
            msc=:black,
            lc=:black,
            grid = false,
            colorbar = false,
            size=(1000, 1000)
        )
    else
        scatter3d(
            eachrow(can.X[:, 1:spatial_downsampling:end])..., 
            marker_z=s,
            # xlabel = "θ₁", ylabel = "θ₂", zlabel = "θ₃",
            xlim=(-1.1, 1.1), ylim=(-1.1, 1.1), zlim=(-1.1, 1.1),
            msc=:black, msa=0.0, msw=0, ms = 8,
            label=nothing, colorbar=false,
            showaxis = false,
            axis=nothing,
            size=(1000, 1000)
        )
    end


    # plot neurons activation in ISO space
    p3 = scatter(
        eachrow(M[:, 1:25:end])...,
        color=:black, alpha=.5, ms=3, 
        showaxis = false,
        axis=nothing,
        xlim = [-p3_extent, p3_extent], ylim = [-p3_extent, p3_extent], zlim = [-p3_extent, p3_extent],
        label=nothing,
        size=(1000, 1000)
    )

    s_embed = S_embedd[
        :, 
        max(hist_frame-100, 1):hist_frame,
    ]
    scatter!(
        eachrow(s_embed)...,
        marker_z = 1:size(s_embed, 2),
        alpha = 1:size(s_embed, 2),
        msa=0, msw=0, colorbar = false,
        label=nothing
    )

    plt = plot(p1, p2, 
        p3, layout=(1,3), 
        # layout=(1,2),
        size=(2200, 800),
        xtickfontsize=16,
        ytickfontsize=16,
        ztickfontsize=16,
        xguidefontsize=16,
        yguidefontsize=16,
        zguidefontsize=16,
        legendfontsize=16,
        right_margin = 12Plots.mm,
        left_margin = 12Plots.mm,
        top_margin = 12Plots.mm,
        bottom_margin = 12Plots.mm,
    )
end