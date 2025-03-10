"""
Load random trajectory simulations for each network and 
measure error and plot it. 
"""

import GeneralAttractors: moving_average
using StatsPlots

# TODO make this work for non standard cover maps

include("settings.jl")
move_to_datadir(supervisor, "PI2")
tag = "PI_random_trajectories"

"""
Measure error as the distance between the true and decoded trajectory
as given by the CAN's coverspace metric.
"""
function measure_error(distance_fn, X, X̄)
    # at each time point, measure the distance between the true and decoded trajectory
    e =  colwise(distance_fn, X', X̄')
    # e[isnan.(e)] .= 0.0

    # at each time point, measure the length of the true trajectory
    l = map(
        t -> distance_fn(X[t, :]', X[t+1, :]), 1:size(X, 1)-1
    ) |> cumsum
    # return moving_average(e, 21)
    return l, e
end


function error_at_checkpoint(e, l, checkpoint)
    # find the first time point where the length of the true trajectory
    # is greater than the checkpoint
    i = findfirst(l .>= checkpoint)

    return if isnothing(i)
        nothing
    else
        isnan(e[i]) ? nothing : (e[i]/checkpoint * 100) # /l * 100  # error in percentage of distance
    end
end




# ---------------------------------------------------------------------------- #
#                                      run                                     #
# ---------------------------------------------------------------------------- #

funky = false
η = 1.5
cover_manifold = "default"

ρ(x) = x

plt = plot(grid=false)

distance_thresholds = Dict(
    "line" => Dict(:low => 1, :medium=>3, :high=>2.5),
    "ring" => Dict(:low => .5, :medium=>1.5, :high=>1.0),
    "torus" => Dict(:low => 1, :medium=>5, :high=>10),
    "sphere" => Dict(:low => 0.5, :medium=>1.5, :high=>2),
    "mobius" => Dict(:low => 1, :medium=>2, :high=>2.0),
    "cylinder" => Dict(:low => 1, :medium=>5, :high=>10),
    "plane" => Dict(:low => 1, :medium=>5, :high=>10),
)

xδ = Dict(
    :low => -1, :medium=>0.0, :high=>1
)

xticksmarks, xticks = [], []
# for (i, network) in enumerate(networks)

network = "torus"
for (i, η) in enumerate((0.0, 1.5, 3.0, 5.0,))
    η > 0 && network != "torus" && continue
    funky == true && network ∉ ("torus", "sphere") && continue
    cover_manifold != "default" && network ∉ ("mobius", ) && continue

    network ∉ ("torus", ) && continue

    # --------------------------------- get data --------------------------------- #
    println(hLine("$(network) - η=$(η)"; style="bold blue"))
    can = make_path_int_can(network)[1]

    filters = Dict(
        :can => network,
        :funky => funky, 
        :noise => η,
        :tag => tag,
        # :cover_manifold => cover_manifold,
    )


    # load data
    meta, data = ProjectSupervisor.fetch(supervisor; filters...)

    # split the different things saved
    trajectories = getfield.(get.(data, "trajectory", nothing), :X)
    decoded_trajectories = get.(data, "decoded", nothing)
    distance_fn = metrics[typeof(can.C.M)]
    @info "Loaded $(length(trajectories)) trajectories"

    errors = map(x -> measure_error(distance_fn, x[1], x[2]), zip(trajectories, decoded_trajectories))
    errors_at_checkpoints = Dict(
        k => filter(
            !isnothing,
            map(e -> error_at_checkpoint(e[2], e[1], distance_thresholds[network][k]), errors) 
        ) 
        for k in keys(distance_thresholds[network])
    )

    # ----------------------------------- plot ----------------------------------- #
    for (j, (k, e)) in enumerate(errors_at_checkpoints)
        k != :high && continue
        @info "Network $network got $(length(e)) errors at checkpoint $k"

        # select 40 entries at random
        e = length(e) > 40 ? e[1:40] : e

        length(e) == 0 && continue
        δ = xδ[k]
        X = zeros(length(e)) .+ 5i .+ δ

        push!(xticksmarks, 5i + δ)
        # push!(xticks, distance_thresholds[network][k])
        push!(xticks, network * " η = $(η)")

        boxplot!(
            X,
            e, 
            label=nothing,
            color = networks_colors[i],
            xticks=(xticksmarks, xticks)
        )
    end

    display(plt)   
    # break
end

display(plt)
save_plot(supervisor, plt, "f5_B_PI_random_trajectories_funky_$(funky)_noise_$(η)_cover_$(cover_manifold)")