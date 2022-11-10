using GeneralAttractors
using GeneralAttractors.ManifoldUtils
using GeneralAttractors.Simulations
import GeneralAttractors.ManifoldUtils: sphere_embedding as φ

using Plots
import LinearAlgebra: norm
import ForwardDiff: jacobian

import GeneralAttractors: by_column

T = 600
traj = Trajectory(S²; T=T)


plt = plot(xlim=[-4, 4], ylim=[-2, 2])


plot!(traj.X[:, 1], traj.X[:, 2])

for i in 1:2:T
    plot!(
        [traj.X[i, 1], traj.X[i, 1]+traj.V[i, 1]],
        [traj.X[i, 2], traj.X[i, 2]+traj.V[i, 2]],
        label=nothing, color=:black
    )
end

plt