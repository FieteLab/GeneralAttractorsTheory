using GeneralAttractors
using Plots
using Term: STACKTRACE_HIDDEN_MODULES
import GeneralAttractors: by_column

STACKTRACE_HIDDEN_MODULES[] = ["Plots"]

include("../networks/mobius.jl")

X = mobiuscan.X
M = by_column(mobius_embedding, X)

idx = 9
plts = []
for i = [1, 3, 6]
    conn = mobiuscan.Ws[i][idx, :] #  .- spherecan_noffset.Ws[i][idx, :]
    p = Plots.scatter3d(
        eachrow(M)...,
        marker_z = conn,
        markersize = 4,
        msa = 0,
        msw = 0,
        alpha = 1,
        label = nothing,
        camera = (90, 0),
        color = :bwr,
        xlim = [0, 1.5],
    )

    scatter3d!(
        [M[1, idx]],
        [M[2, idx]],
        [M[3, idx]],
        ms = 8,
        color = :black,
    )

    push!(plts, p)
end

plot(plts..., layout = (3, 1), size = (800, 1000), colorbar = nothing)
