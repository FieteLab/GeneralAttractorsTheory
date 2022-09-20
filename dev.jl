using GeneralAttractors
using Plots

import GeneralAttractors.Simulations: RandomChunk

RandomChunk(simulation; duration=10000, μ₀=1.0, σ=3) |> plot