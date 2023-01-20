"""
Definition of parameters and settins used throughout the analysis/plots.
"""

using Plots, DataFrames, Term, Statistics
gr()
# pyplot()
install_term_stacktrace(; hide_frames=false)

using GeneralAttractors.Simulations
using GeneralAttractors.Analysis
using GeneralAttractors.Kernels
import GeneralAttractors: torus_maker, sphere_maker, mobius_maker, cylinder_maker, plane_maker
import GeneralAttractors: torus_embedding, identity_embedding, mobius_embedding, cylinder_embedding, plane_embedding
using GeneralAttractors.ProjectSupervisor
import GeneralAttractors: by_column

datadir = "/Users/federicoclaudi/Desktop/GeneralAttractors/data"

supervisor = Supervisor("GeneralAttractorsTheory")
set_datadir(supervisor, datadir)

# ---------------------------------- general --------------------------------- #
τ = 5.0
b₀ = 1.0
dt = 0.5

networks = ("torus", "sphere", "mobius", "cylinder", "plane")

network_makers = Dict(
    "torus" => torus_maker,
    "sphere" => sphere_maker,
    "mobius" => mobius_maker,
    "cylinder" => cylinder_maker,
    "plane" => plane_maker,
)

embeddings = Dict(
    "torus" => torus_embedding,
    "sphere" => identity_embedding,  # because the neurons are already "embedded" on the sphere
    "mobius" => mobius_embedding,
    "cylinder" => cylinder_embedding,
    "plane" => plane_embedding,
)


# ---------------------------- ANALYSIS paramters ---------------------------- #
dimred_3d_params = AnalysisParameters(
    max_nPC = 50,
    pca_pratio = 0.9999,
    n_isomap_dimensions = 3,
    isomap_k = 10,
    isomap_downsample = 50,  # time downsamplin
)

dimred_10d_params = AnalysisParameters(
    max_nPC = 50,
    pca_pratio = 0.9999,
    n_isomap_dimensions = 10,
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
            :σ => 10:δ:100,
        ),
        :DoE => Dict(
            :a => 4:δ:20,
            :λ => 20:δ:100,
        ),
        :local_global => Dict(
            :α => 0.8:δ:5,
            :σ => 150:δ:200,
        ),
        :constant => Dict(
            :σ => 1.5:δ:3.0,
            :β_minus => -1:δ:-.1
        ),
    ),

"mobius" => Dict(
        :mexican_hat => Dict(
            :α => 3:δ:16,
            :σ => 2:δ:30,
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
            :σ => 2:δ:30,
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
            :σ => 2:δ:30,
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
            :σ => 2:δ:30,
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