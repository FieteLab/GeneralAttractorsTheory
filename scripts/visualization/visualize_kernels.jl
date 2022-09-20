using GeneralAttractors
using Plots

using GeneralAttractors.Kernels

plot(
    plot(MexicanHatKernel(), title="Mexican hat"),
    plot(DiffOfExpKernel(), title="Difference of exp"),
)