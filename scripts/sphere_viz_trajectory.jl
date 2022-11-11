using GeneralAttractors
using GeneralAttractors.ManifoldUtils
using GeneralAttractors.Simulations
import GeneralAttractors.ManifoldUtils: sphere_embedding as φ

using Plots
import LinearAlgebra: norm
import ForwardDiff: jacobian

import GeneralAttractors: by_column

T = 20000
traj = Trajectory(S²; T=T, σ=[0.0, 0.0, 0.0])


plt = plot(xlim=[-1.1, 1.1], ylim=[-1.1, 1.1], zlim=[-1.1, 1.1])


plot!(eachcol(traj.X)...)

# for i in 1:2:T
#     plot!(
#         [traj.X[i, 1], traj.X[i, 1]+traj.V[i, 1]],
#         [traj.X[i, 2], traj.X[i, 2]+traj.V[i, 2]],
#         label=nothing, color=:black
#     )
# end

plt