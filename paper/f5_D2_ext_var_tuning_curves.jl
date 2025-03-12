# using StatsPlots
# using Statistics

# include("settings.jl")
# move_to_datadir(supervisor, "TuningCurves")
tag = "PI_systematic_trajectories"



# ---------------------------------------------------------------------------- #
#                                Load the data                                 #
# ---------------------------------------------------------------------------- #

function load_trajectory_data(network, η=1.5; funky=false)
    # Set up filters for data loading
    filters = Dict(
        :can => network,
        :funky => funky, 
        :noise => η,
        :tag => tag,
    )

    # Load data
    meta, data = ProjectSupervisor.fetch(supervisor; filters...)

    # Extract trajectories and neural activities
    trajectories = getfield.(get.(data, "h", nothing), :x_M)  # n_dim x n_steps
    activities = getfield.(get.(data, "h", nothing), :S)  # n_neurons x n_cans x n_steps

    return trajectories, activities
end

# ---------------------------------------------------------------------------- #
#                                    Plotting                                  #
# ---------------------------------------------------------------------------- #

function plot_neuron_tuning(trajectories, activities, network; neuron_idx=1, n_bins=51)
    # Concatenate all trajectories and activities
    all_positions = hcat(trajectories...)  # n_dim x (n_steps * n_trajectories)
    all_activities = vcat([act[neuron_idx, 2, :] for act in activities]...)
    
    # Extract x and y coordinates
    x_coords = all_positions[1, :]

    if network in ("ring", "line")
        # Create a scatter plot for 1D tuning curves
        plt = scatter(
            x_coords, 
            all_activities, 
            marker_z=all_activities,
            xlabel="x[1]", 
            ylabel="Activity", 
            title="Neuron $neuron_idx Tuning Curve",
            markersize=5,
            alpha=0.5
        )
    else
        y_coords = all_positions[2, :]

        # Create bins for x and y coordinates
        x_edges = range(minimum(x_coords), maximum(x_coords), n_bins)
        y_edges = range(minimum(y_coords), maximum(y_coords), n_bins)
        
        # Initialize matrix to store binned activities
        activity_grid = zeros(n_bins-1, n_bins-1)
        count_grid = zeros(n_bins-1, n_bins-1)
        
        # Bin the activities
        for (x, y, a) in zip(x_coords, y_coords, all_activities)
            # add noise to the coordinates
            x += randn() * 0.1
            y += randn() * 0.1

            i = searchsortedfirst(x_edges, x) - 1
            j = searchsortedfirst(y_edges, y) - 1
            
            # Ensure indices are within bounds
            if 1 <= i < n_bins && 1 <= j < n_bins
                activity_grid[i, j] += a
                count_grid[i, j] += 1
            end
        end
        
        # Average the activities in each bin
        activity_grid ./= count_grid .+ eps()  # add eps to avoid division by zero
        
        # Create heatmap
        plt = heatmap(
            x_edges[1:end-1], 
            y_edges[1:end-1], 
            activity_grid',  # transpose to match coordinate system
            color=:viridis,
            title="Neuron $neuron_idx Tuning Curve",
            xlabel="x position",
            ylabel="y position",
            colorbar_title="Mean Activity",
            aspect_ratio=:equal
        )
    end
    
    return plt
end

# ---------------------------------------------------------------------------- #
#                                      Run                                     #
# ---------------------------------------------------------------------------- #

# Load data for the desired network
network = "ring"  # or "line", "plane", etc.
η = 0.0
trajectories, activities = load_trajectory_data(network, η)

# Print some debug info about the loaded data
@info "Number of trajectories:" length(trajectories)
@info "First trajectory size:" size(first(trajectories))
@info "First activity tensor size:" size(first(activities))

# Create a figure with multiple neurons
n_neurons = 4
for (k, i) in enumerate(1:n_neurons)
    neuron_idx = rand(1:size(activities[1], 1))
    p = plot_neuron_tuning(trajectories, activities, network; neuron_idx=neuron_idx)
    display(p)
    save_plot(supervisor, p, "f5_D_ext_var_tuning_curves_$(network)_$(i)")
end

