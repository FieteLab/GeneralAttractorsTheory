using GeneralAttractors
using GeneralAttractors.ManifoldUtils
using GeneralAttractors.Simulations
import GeneralAttractors.ManifoldUtils: sphere_embedding as φ

using Plots
import LinearAlgebra: norm
import ForwardDiff: jacobian

import GeneralAttractors: by_column

T = 2500
x₀ = [0, 1, 0]
traj = Trajectory(S²; T = T, σ = [0.0, 0.0, 0.0], x₀ = x₀, scale = 0.05, modality = :rand)


plt = Plots.plot(xlim = [-1.1, 1.1], ylim = [-1.1, 1.1], zlim = [-1.1, 1.1])


Plots.plot3d!(
    [-1, 1],
    [-1, -1],
    [0, 0],
    lw = 2,
    label = nothing,
    color = :black,
    alpha = 0.2,
)
Plots.plot3d!(
    [-1, -1],
    [-1, 1],
    [0, 0],
    lw = 2,
    label = nothing,
    color = :black,
    alpha = 0.2,
)
Plots.plot3d!([-1, 1], [1, 1], [0, 0], lw = 2, label = nothing, color = :black, alpha = 0.2)
Plots.plot3d!([1, 1], [-1, 1], [0, 0], lw = 2, label = nothing, color = :black, alpha = 0.2)

Plots.plot!(eachcol(traj.X)..., lw = 3, color = :black, label = "traj")

# for i in 1:2:T
#     plot!(
#         [traj.X[i, 1], traj.X[i, 1]+traj.V[i, 1]],
#         [traj.X[i, 2], traj.X[i, 2]+traj.V[i, 2]],
#         label=nothing, color=:black
#     )
# end

plt
