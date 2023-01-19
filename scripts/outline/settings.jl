"""
Definition of parameters and settins used throughout the analysis/plots.
"""

using Plots
# gr()
pyplot()

using GeneralAttractors.Simulations
using GeneralAttractors.Analysis
using GeneralAttractors.Kernels
import GeneralAttractors: torus_maker
import GeneralAttractors: torus_embedding, sphere_embedding, mobius_embedding
using GeneralAttractors.ProjectSupervisor
import GeneralAttractors: by_column

datadir = "/Users/federicoclaudi/Dropbox (Personal)/Postdoc/GeneralAttractors/data"

supervisor = Supervisor("GeneralAttractorsTheory")
set_datadir(supervisor, datadir)

# ---------------------------------- general --------------------------------- #
τ = 5.0
b₀ = 1.0


network_makers = Dict(
    "torus" => torus_maker,
#     "sphere" => sphere_maker,
#     "mobius" => mobius_maker,
)
embeddings = Dict(
    "torus" => torus_embedding,
    "sphere" => sphere_embedding,
    "mobius" => mobius_embedding,
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
            :α => 0.3:δ:2.5,
            :σ => 5:δ:100,
        ),
        :DoE => Dict(
            :a => 1.0:δ:2.0,
            :λ => 5:δ:80,
        ),
        :local_global => Dict(
            :α => .1:δ:3,
            :σ => 150:δ:300,
        ),
        :constant => Dict(
            :σ => 1.5:δ:3.0,
            :β_minus => -1:δ:-.1
        ),
    )
)




include("utils.jl")