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
import GeneralAttractors: torus_maker, sphere_maker, mobius_maker, cylinder_maker
import GeneralAttractors: torus_embedding, sphere_embedding, mobius_embedding, cylinder_embedding
using GeneralAttractors.ProjectSupervisor
import GeneralAttractors: by_column

datadir = "/Users/federicoclaudi/Desktop/GeneralAttractors/data"

supervisor = Supervisor("GeneralAttractorsTheory")
set_datadir(supervisor, datadir)

# ---------------------------------- general --------------------------------- #
τ = 5.0
b₀ = 1.0
dt = 0.5


network_makers = Dict(
    "torus" => torus_maker,
    "sphere" => sphere_maker,
    "mobius" => mobius_maker,
    "cylinder" => cylinder_maker,
)

embeddings = Dict(
    "torus" => torus_embedding,
    "sphere" => sphere_embedding,
    "mobius" => mobius_embedding,
    "cylinder" => cylinder_embedding,
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
    )
)




include("utils.jl")