using Plots
using GeneralAttractors
using GeneralAttractors.ManifoldUtils
import GeneralAttractors: by_column

include("../networks/mobius.jl")

"""
Plot the Mobius lattice in 3d colored by the value
of each manifold coordinate.
"""

X = mobiuscan.X

M = by_column(mobius_embedding, X)

x1, x2 = X[1, :], X[2, :]



plot(
    scatter3d(eachrow(M)..., marker_z = x1, label = nothing, title = "X"),
    scatter3d(eachrow(M)..., marker_z = x2, label = nothing, title = "Y"),
    layout = (2, 1),
    size = (400, 800),
    camera = (40, 30),
    msw = 0.25,
)
