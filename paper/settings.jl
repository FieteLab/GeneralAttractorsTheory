"""
Definition of parameters and settins used throughout the analysis/plots.
"""

using Plots, DataFrames, Term, Statistics, LinearAlgebra, ObjectivePaths, Distances
import MyterialColors: Palette, green_dark, deep_purple, indigo, salmon, salmon_dark

gr()   # for fast plotting
# pyplot()  # for better plots quality

install_term_stacktrace(; hide_frames=true)

using GeneralAttractors.Simulations
using GeneralAttractors.Analysis
using GeneralAttractors.Kernels
import GeneralAttractors: torus_maker, sphere_maker, mobius_maker, cylinder_maker, plane_maker
import GeneralAttractors: torus_embedding, identity_embedding, mobius_embedding, cylinder_embedding, plane_embedding, sphere_embedding
using GeneralAttractors.ProjectSupervisor
import GeneralAttractors: by_column, MobiusEuclidean, SphericalDistance, moving_average
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: remove_jumps_from_trajectory

datadir = "/Users/federicoclaudi/Desktop/GeneralAttractors/data"

supervisor = Supervisor("GeneralAttractorsTheory")
set_datadir(supervisor, datadir)

# ---------------------------------- general --------------------------------- #
τ = 5.0
b₀ = 1.0
dt = 0.5

networks = (
    "ring",
    "torus", 
    "sphere", 
    "mobius", 
    "cylinder",
    "plane"
)

cans = Dict(
    "ring" => ring_maker(:default),
    "torus" => torus_maker(:default), 
    "sphere" => sphere_maker(:default),
    "mobius" => mobius_maker(:default),
    "cylinder" => cylinder_maker(:default),
    "plane" => plane_maker(:default),
)

network_makers = Dict(
    "ring" => ring_maker,
    "torus" => torus_maker,
    "sphere" => sphere_maker,
    "mobius" => mobius_maker,
    "cylinder" => cylinder_maker,
    "plane" => plane_maker,
)

embeddings = Dict(
    "ring" => ring_embedding,
    "torus" => torus_embedding,
    "sphere" => identity_embedding,  # because the neurons are already "embedded" on the sphere
    "mobius" => mobius_embedding,
    "cylinder" => cylinder_embedding,
    "plane" => plane_embedding,
)

metrics = Dict(
    Ring => PeriodicEuclidean([2π]),
    Manifoldℝ² => Euclidean(),
    Cylinder => PeriodicEuclidean([2π, Inf]),
    Torus => PeriodicEuclidean([2π, 2π]),
    Sphere => SphericalDistance(),
    Mobius => MobiusEuclidean(),
)


max_path_int_vel = Dict(
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
    :grid => false,
    :legendfontsize=>16,
    :right_margin => 12Plots.mm,
    :left_margin => 12Plots.mm,
    :top_margin => 12Plots.mm,
    :bottom_margin => 12Plots.mm,
    :dpi => 300,
)

# ---------------------------- ANALYSIS paramters ---------------------------- #
dimred_3d_params = AnalysisParameters(
    max_nPC = 400,  # max num of PCs
    pca_pratio = 0.999999999999999999,       # fraction of variance explained
    n_isomap_dimensions = 3,
    isomap_k = 10,
    isomap_downsample = 50,
)

dimred_10d_params = AnalysisParameters(
    max_nPC = 200,
    pca_pratio = 0.9999999,
    n_isomap_dimensions = 10,
    isomap_k = 10,
    isomap_downsample = 50,  # time downsamplin
)


dimred_50d_params = AnalysisParameters(
    max_nPC = 200,
    pca_pratio = 0.9999999,
    n_isomap_dimensions = 50,
    isomap_k = 10,
    isomap_downsample = 50,  # time downsamplin
)




# --------------------------- kernel values ranges --------------------------- #

"""
For each kernel type and CAN architecture 
"""

kernels = (;
   :mexican_hat => MexicanHatKernel, 
   :DoE => DiffOfExpKernel,
   :local_global => LocalGlobalKernel,
   :constant => ConstantKernel,
)


δ = 0.01
kernels_parameters_range = Dict(
    "torus" => Dict(
        :mexican_hat => Dict(
            :α => 3:δ:16,
            :σ => 3:δ:80,
        ),
        :DoE => Dict(
            :a => 4:δ:20,
            :λ => 10:δ:100,
        ),
        :local_global => Dict(
            :α => 0.8:δ:8,
            :σ => 100:δ:200,
        ),
        :constant => Dict(
            :σ => 1.0:δ:2.0,
            :β_minus => -1:δ:-.1
        ),
    ),

"mobius" => Dict(
        :mexican_hat => Dict(
            :α => 3:δ:16,
            :σ => 1:δ:30,
        ),
        :DoE => Dict(
            :a => 4:δ:20,
            :λ => 4:δ:30,
        ),
        :local_global => Dict(
            :α => 0.8:δ:8,
            :σ => 50:δ:75,
        ),
        :constant => Dict(
            :σ => 0.5:δ:1.5,
            :β_minus => -1:δ:-.1
        ),
    ),

"cylinder" => Dict(
        :mexican_hat => Dict(
            :α => 3:δ:16,
            :σ => 1:δ:30,
        ),
        :DoE => Dict(
            :a => 4:δ:20,
            :λ => 4:δ:30,
        ),
        :local_global => Dict(
            :α => 0.8:δ:8,
            :σ => 50:δ:75,
        ),
        :constant => Dict(
            :σ => 0.5:δ:1.5,
            :β_minus => -1:δ:-.1
        ),
    ),

    "plane" => Dict(
        :mexican_hat => Dict(
            :α => 3:δ:16,
            :σ => 1:δ:30,
        ),
        :DoE => Dict(
            :a => 4:δ:20,
            :λ => 4:δ:30,
        ),
        :local_global => Dict(
            :α => 0.8:δ:8,
            :σ => 50:δ:75,
        ),
        :constant => Dict(
            :σ => 0.5:δ:1.0,
            :β_minus => -1:δ:-.1
        ),
    ),

    "sphere" => Dict(
        :mexican_hat => Dict(
            :α => 3:δ:16,
            :σ => 1:δ:30,
        ),
        :DoE => Dict(
            :a => 4:δ:20,
            :λ => 4:δ:30,
        ),
        :local_global => Dict(
            :α => 2:δ:10,
            :σ => 50:δ:75,
        ),
        :constant => Dict(
            :σ => 0.8:δ:1.2,
            :β_minus => -1:δ:-.1
        ),
    ),
)




include("utils.jl")