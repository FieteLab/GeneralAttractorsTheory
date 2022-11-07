using GeneralAttractors
using Plots

using GeneralAttractors.Kernels

plot(
    plot(MexicanHatKernel(α=0.5,  σ=1.0), title="Mexican hat"),
    plot(DiffOfExpKernel(), title="Difference of exp"),
    plot(LocalGlobalKernel(), title="Local ex global inh"),
)