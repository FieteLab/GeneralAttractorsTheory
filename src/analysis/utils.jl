function ∑(x)
    n, _, m = size(x)
    real.(reshape(mean(x, dims = 2), (n, m)))
end

"""
Get the average on mfld bump velocity for a simultion
"""
function get_bump_speed(can, fld_name::String, sim_name::String)::Float64
    # load state history
    history = load_simulation_history(fld_name, sim_name)
    # @info "h" history history.S history.Δt

    s = ∑(history.S)[:, 10:end]

    # get peak location speed
    peak_location = hcat(
        map(
            st -> decode_peak_location(st, can), 
            eachcol(s)
            )...
        )

    on_mfld_speed = map(
        i -> can.metric(
            peak_location[:, i], peak_location[:, i-1]
            ), 
        2:size(peak_location, 2)
    )
    average_speed = sum(on_mfld_speed)/(length(on_mfld_speed)*history.Δt)  # tot displacement over time
    return average_speed 
end