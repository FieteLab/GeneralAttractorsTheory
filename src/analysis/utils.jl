function ∑(x)
    n, _, m = size(x)
    real.(reshape(sum(x, dims = 2), (n, m)))
end

"""
Get the average on mfld bump velocity for a simultion
"""
function get_bump_speed(can, fld_name::String, sim_name::String)
    # load state history
    history = load_simulation_history(fld_name, sim_name)
    return get_bump_speed(can,history)
end


function get_bump_speed(can, history)
    s = ∑(history.S)[:, 20:end]

    # get peak location speed
    peak_location = hcat(map(st -> decode_peak_location(st, can), eachcol(s))...)

    on_mfld_speed = map(
        i -> i==1 ? 0 : can.metric(peak_location[:, i], peak_location[:, i-1]),
        1:size(peak_location, 2),
    ) / history.Δt
    average_speed = sum(on_mfld_speed) / length(on_mfld_speed) 
    return on_mfld_speed, average_speed
end
