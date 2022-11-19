using Plots


using GeneralAttractors
using GeneralAttractors.Simulations
using Interpolations
using Term
using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp, by_column
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: plot_trajectory_and_decoded
import GeneralAttractors.Simulations: decode_peak_location
using GeneralAttractors.Analysis.ManifoldAnalysis: population_average
import GeneralAttractors: torus_embedding, sphere_embedding, mobius_embedding
using ManifoldsBase, Manifolds

# include("../scripts/networks/torus.jl")
# include("../scripts/networks/sphere.jl")
# include("../scripts/networks/mobius.jl")


sim_fld = "sphere"
sim_name = "decoding"

if sim_fld == "torus"
    can = toruscan
    φ = torus_embedding
elseif sim_fld == "sphere"
    can = spherecan
    φ = sphere_embedding
else
    can = mobiuscan
    φ = mobius_embedding
end


sim = "$(sim_name)_$(sim_fld)"
history = load_simulation_history(sim_fld, sim*"_history")
S = population_average(history)

X̄ = load_data(sim_fld, sim_name*"_"*sim_fld*"_sim_decoded_X")
X = load_data(sim_fld, sim_name*"_"*sim_fld*"_sim_trajectory_X")

# ------------- plot bump trajectory on the maniold lattice in 3D ------------ #
# TODO this doesn't actually show the bump but the bump through the decoder to fix in the future
mfld_position = by_column(c -> decode_peak_location(c, can), S)
t = 1:size(mfld_position, 2)
fm = linear_interpolation(t, mfld_position[1, :])
fn = linear_interpolation(t, mfld_position[2, :])

T = 1:.01:size(mfld_position, 2)
m = fm.(T)
n = fn.(T)

mfld_position = Matrix(hcat(m, n)')


embed_position = by_column(φ, mfld_position)


# scatter3d(eachrow(embed_position)..., ms=8, alpha=.2, color="#BABDD3",
# plot3d(eachrow(embed_position)..., lw=5, color="#9CA1C7",
# msa=0, msw=0,
# grid=false,
# xticks=nothing,
# yticks=nothing,
# zticks=nothing,
# xlim=[-1.5, 1.5],
# ylim=[-1.5, 1.5],
# zlim=[-1.5, 1.5],
# camera=(45, 55)
# )


# -------------------- plot decoded and actual trajectory -------------------- #
plt = plot(
    xlabel="x₁",
    ylabel="x₂",
    grid=false,
    camera=(90, 60)
)

plot!(eachcol(X)..., lw=4, color=:red, label=nothing)
plot!(eachcol(X̄)..., lw=3, color=:black, label=nothing)