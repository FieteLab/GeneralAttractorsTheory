using GeneralAttractors
using Plots

using GeneralAttractors.Kernels

plot(
    plot(MexicanHatKernel(α = 0.01, σ = 0.1, β=0), title = "Mexican hat"),
    # plot(DiffOfExpKernel(), title = "Difference of exp"),
    # plot(LocalGlobalKernel(σ=0.33), title = "Local ex global inh"),
)
