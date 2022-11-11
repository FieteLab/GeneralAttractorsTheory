using GeneralAttractors
using Plots

using GeneralAttractors.Kernels

plot(
    plot(MexicanHatKernel(α = 0.5, σ = 1.0, β=1), title = "Mexican hat"),
    plot(DiffOfExpKernel(), title = "Difference of exp"),
    plot(LocalGlobalKernel(σ=0.33), title = "Local ex global inh"),
)
