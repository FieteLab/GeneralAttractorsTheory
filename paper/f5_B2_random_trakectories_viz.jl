"""
Load random trajectory simulations for each network and 
measure error and plot it. 
"""

import GeneralAttractors: moving_average
using StatsPlots

# TODO make this work for non standard cover maps

include("settings.jl")
move_to_datadir(supervisor, "PI")
tag = "PI_random_trajectories"

"""
Measure error as the distance between the true and decoded trajectory
as given by the CAN's coverspace metric.
"""
function measure_error(distance_fn, X, X̄)
    e =  colwise(distance_fn, X', X̄')
    
    # replace nans with 0
    e[isnan.(e)] .= 0.0
    return moving_average(e, 21)
end




# ---------------------------------------------------------------------------- #
#                                      run                                     #
# ---------------------------------------------------------------------------- #

funky = true
η = 0.0

ρ(x) = x

plt = plot()

for (i, network) in enumerate(networks)

    η > 0 && network != "torus" && continue
    funky == true && network ∉ ("torus", "sphere") && continue

    # network ∉ ("plane", ) && continue

    # --------------------------------- get data --------------------------------- #
    println(hLine("$(network) - η=$(η)"; style="bold blue"))
    can = make_path_int_can(network)[1]

    filters = Dict(
        :can => network,
        :funky => funky, 
        :noise => η,
        :tag => tag,
    )


    # load data
    meta, data = ProjectSupervisor.fetch(supervisor; filters...)

    # split the different things saved
    trajectories = getfield.(get.(data, "trajectory", nothing), :X)
    decoded_trajectories = get.(data, "decoded", nothing)
    distance_fn = metrics[typeof(can.C.M)]

    errors = map(x -> measure_error(distance_fn, x[1], x[2])[end], zip(trajectories, decoded_trajectories))
    errors = hcat(errors...)[1, :]

    # ----------------------------------- plot ----------------------------------- #

    # do a box plot to show the distribution of final error across repeats
    X = zeros(length(errors)) .+ i
    violin!(
        X,
        errors, 
        label=network
    )
    boxplot!(
        X,
        errors, 
        label=nothing
    )
    dotplot!(
        X,
        errors, 
        label=nothing
    )

    
    display(plt)
end

save_plot(supervisor, plt, "f5_B_PI_random_trajectories_funky_$(funky)_noise_$(η)")