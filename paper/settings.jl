"""
Definition of parameters and settins used throughout the analysis/plots.
"""

using Plots, DataFrames, Term, Statistics, LinearAlgebra, ObjectivePaths, Distances
import MyterialColors: Palette, green_dark, deep_purple, indigo, salmon, salmon_dark
import Plots: grid
using ProjectSupervisor
import ProjectSupervisor

gr()   # for fast plotting
# pyplot()  # for better plots quality

install_term_stacktrace(; hide_frames=false)

using GeneralAttractors.Simulations
using GeneralAttractors.Analysis
using GeneralAttractors.Kernels
import GeneralAttractors: torus_maker, sphere_maker, mobius_maker, cylinder_maker, plane_maker, ring_maker, line_maker
import GeneralAttractors: torus_embedding, identity_embedding, mobius_embedding, cylinder_embedding, plane_embedding, sphere_embedding, ring_embedding, line_embedding
import GeneralAttractors: by_column, MobiusEuclidean, SphericalDistance, moving_average
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: remove_jumps_from_trajectory
import GeneralAttractors.Simulations: plot_trajectory_and_decoded

datadir = "/Users/federicoclaudi/Desktop/GeneralAttractors/data"

supervisor = Supervisor("GeneralAttractorsTheory"; 
            plotsdir="/Users/federicoclaudi/Dropbox (Personal)/Postdoc/Writings/GeneralAttractors/plots")
set_datadir(supervisor, datadir)

# ---------------------------------- general --------------------------------- #
τ = 5.0
b₀ = 1.0
dt = 0.5

networks = (
    "line",
    "ring",
    "torus", 
    "sphere", 
    "mobius", 
    "cylinder",
    "plane"
)



network_makers = Dict(
    "line" => line_maker,
    "ring" => ring_maker,
    "torus" => torus_maker,
    "sphere" => sphere_maker,
    "mobius" => mobius_maker,
    "cylinder" => cylinder_maker,
    "plane" => plane_maker,
)

embeddings = Dict(
    "line" => line_embedding,
    "ring" => ring_embedding,
    "torus" => torus_embedding,
    "sphere" => identity_embedding,  # because the neurons are already "embedded" on the sphere
    "mobius" => mobius_embedding,
    "cylinder" => cylinder_embedding,
    "plane" => plane_embedding,
)

metrics = Dict(
    Line => Euclidean(),
    Ring => PeriodicEuclidean([2π]),
    Manifoldℝ² => Euclidean(),
    Cylinder => PeriodicEuclidean([2π, Inf]),
    Torus => PeriodicEuclidean([2π, 2π]),
    Sphere => SphericalDistance(),
    Mobius => MobiusEuclidean(),
)


max_path_int_vel = Dict(
    "line" => 0.04,
    "ring" => 0.04,
    "plane" => 0.04,
    "cylinder" => 0.04,
    "torus" => 0.04,
    "sphere" => 0.004,
    "mobius" => 0.02,
)

# ------------------------------ PLOTTING params ----------------------------- #
networks_colors = getfield.(Palette(green_dark, deep_purple; N = length(networks)).colors, :string)


plot_font_size_kwargs = Dict(
    :xtickfontsize=>16,
    :ytickfontsize=>16,
    :ztickfontsize=>16,
    :xguidefontsize=>16,
    :yguidefontsize=>16,
    :zguidefontsize=>16,
    :legendfontsize=>16,
    :right_margin => 12Plots.mm,
    :left_margin => 12Plots.mm,
    :top_margin => 12Plots.mm,
    :bottom_margin => 12Plots.mm,
    :dpi => 300,
)

plot_remove_axes = Dict(
    :grid => false,
    :showaxis => false,
    :axis=>nothing,
)

# ---------------------------- ANALYSIS paramters ---------------------------- #
dimred_3d_params = AnalysisParameters(
    max_nPC = 400,  
    pca_pratio = 0.999999999999999999,     
    n_isomap_dimensions = 3,
    isomap_k = 10,
    isomap_downsample = 10,
)

dimred_10d_params = AnalysisParameters(
    max_nPC = 200,
    pca_pratio = 0.9999999,
    n_isomap_dimensions = 10,
    isomap_k = 10,
    isomap_downsample = 10,  
)


dimred_50d_params = AnalysisParameters(
    max_nPC = 200,
    pca_pratio = 0.9999999,
    n_isomap_dimensions = 50,
    isomap_k = 10,
    isomap_downsample = 10,  
)




intrinsic_dimensionality_prms = AnalysisParameters(
    intrinsic_d_nseeds = 250,
    intrinsic_d_pratio = 0.75,
    intrinsic_d_neighborhood_size = 500;
)

include("utils.jl")