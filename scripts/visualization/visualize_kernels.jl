using GeneralAttractors
using Plots

using GeneralAttractors.Kernels

plot(
    plot(MexicanHatKernel(α = 0.1, σ = 0.25), title = "Mexican hat"),
    plot(DiffOfExpKernel(), title = "Difference of exp"),
)
