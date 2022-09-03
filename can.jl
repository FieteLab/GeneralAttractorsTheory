"""
Continuous Attractor Network for grid cells. 

Based on Burak & Fiete 2009.

"""

import Base.Iterators: product as ×  # cartesian product
using Plots
using Term.Progress
import Term

# parameters
n = 64  # √ of number of neurons
f(x) = max(0, x)  # relu
τ = 10  # ms
dt = 0.5 # ms

a = 1.0
α = 0.10315
l = 2.0
λ = 13
β = 3/λ^2
γ = 1.05 * β


""" compute weight for an entry in the connectivity matrix """
W₀(x::Float64) = a * exp(-γ*abs(x)^2) - exp(-β*abs(x)^2)
W₀(x::Vector{Float64}) = W₀.(x) |> sum



function make_network()
    # specify the orientation of each neuron (2x2 tiled square of orientations)
    orientations = [[1, 2] [3, 4]] # 1=North, 2=East, 3=Sout, 4=West
    Θ = repeat(orientations, Int(n/2), Int(n/2))  # n × n mtx
    Θ = reshape(Θ, n^2)  # n^2 vector
    

    # get the position of every neuron in the lattice
    n̂ = Int(n/2)
    X::Matrix{Vector} = [[x...] for x in ((-n̂:(n̂-1)) × (-n̂:(n̂-1)))]

    # make recurrent weight matrix
    E = Dict{Int, Vector}(1=>[0.0, 1.0], 2=>[1.0, 0.0], 3=>[0.0, -1.0], 4=>[-1.0, 0.0])
    W::Matrix{Float64} = zeros(Float64, n^2, n^2)
    
    pbar = ProgressBar(columns=:detailed, expand=false)
    Term.Progress.with(pbar) do
        ijob = addjob!(pbar; description="Computing W", N=n^2)
        for i in 1:n^2
            for j in 1:n^2   
                xᵢ, xⱼ = X[i], X[j]
                ê = E[Θ[i]]
                x = evaluate(metric, xᵢ - (l .* ê), xⱼ)
                W[i, j] = W₀(x)
            end
            update!(ijob)
        end
    end

    return Θ, W, X
end


Θ, W, X = make_network()


plots = []
for i in (1, 15, 2200, 3499)
    weights = reshape(W[i, :], n, n)
    p = heatmap(weights,  colorbar=nothing, aspect_ratio=:equal, xlabel="neuron #", ylabel="neruon #", title= "neuron $i")

    x, y = X[i] .+ n/2 .+ 1
    scatter!([y], [x],  ms=10, color="green", label=nothing)
    push!(plots, p)
end
plot(plots..., size=(500, 500))
