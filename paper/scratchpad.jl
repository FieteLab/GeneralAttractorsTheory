include("settings.jl")


can = network_makers["ring"](:single)


plot(can.W[100, :], 
    fillrange = -1,
    fillalpha = 0.1,
    fillcolor = :red
) |> display