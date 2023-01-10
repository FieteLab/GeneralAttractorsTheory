using Plots
using GeneralAttractors

include("../networks/ring.jl")

plots = []
for i in (100, 200)
    p = plot(ringcan.Ws[1][i, :], lw = 2, color = :red)
    plot!(ringcan.Ws[2][i, :], lw = 2, color = :green)
    push!(plots, p)
end

plot(plots...)
