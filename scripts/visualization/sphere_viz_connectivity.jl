using GeneralAttractors
using Plots


include("../networks/sphere.jl")

spherecan_noffset =
    CAN("sphere", cover, (m, m), I, X, d_s, k_s; offset_size = 0, offsets = offsets, Ω = Ω)

idx = 1000
W = spherecan.Ws[1][idx, :]

# plt = GLMakie.scatter(eachrow(spherecan.X)..., markersize=50, color=W, alpha=0.1)

# for i in 1:100:size(spherecan.X, 2)
#     x = spherecan.X[:, i]

#     for (j, color) in zip((1, 3, 6),(:red, :green, :blue))
#         v = spherecan.offsets[j].ψ(x) .* 0.1
#         arrows!([[y] for y in x]..., [[y] for y in v]..., color=color)
#     end
# end

# println("done")
# plt

plts = [plot(k_s; σ = 0.05)]
for i = 1:2:5
    conn = spherecan.Ws[i][idx, :] #  .- spherecan_noffset.Ws[i][idx, :]
    p = Plots.scatter3d(
        eachrow(spherecan.X)...,
        marker_z = conn,
        clims = (-1, 0),
        markersize = 4,
        msa = 0,
        msw = 0,
        alpha = 1,
        label = nothing,
        camera = (90, 15),
        color = :bwr,
        xlim = [0, 1.5],
    )

    scatter3d!(
        [spherecan.X[1, idx]],
        [spherecan.X[2, idx]],
        [spherecan.X[3, idx]],
        ms = 8,
        color = :black,
    )

    push!(plts, p)
end

plot(plts..., layout = (1, 4), size = (1000, 500), colorbar = nothing)
