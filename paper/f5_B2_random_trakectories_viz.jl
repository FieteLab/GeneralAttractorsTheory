"""
Load random trajectory simulations for each network and 
measure error and plot it. 
"""

# TODO make this work for non standard cover maps

include("settings.jl")
move_to_datadir(supervisor, "PI")


"""
Measure error as the distance between the true and decoded trajectory
as given by the CAN's coverspace metric.
"""
function measure_error(can, X, X̄)
    d = metrics[typeof(can.C.M)]

    e =  [max(d(X[i, :], X̄[i, :]), eps()) for i in 1:size(X, 1)]
    # replace nans with 0
    e[isnan.(e)] .= 0.0
    return e
end




# ------------------------------------ run ----------------------------------- #
network = "torus"
funky = false
η = 0.0

filters = Dict(
    :can => network,
    :funky => funky, 
    :η => η,
)


# load data
meta, data = ProjectSupervisor.fetch(supervisor; filters...)
print(first(meta, 5))

# TODO for each run, load the decoded and the actual trajectory 
# TODO compute error for each run and plot mean/std across runs
