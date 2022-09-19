using Flux
using Term
using Term.Progress
using FluxTraining
using Statistics
using Plots
using Parameters

import Base.Iterators: product as ×

using GeneralAttractors
include("manifolds.jl")

"""
This includes some tests on using small auto encoder 
networks to create global curvliniear coordinate systems
on manifolds embedded in higher dimensional spaces.
"""


# ----------------------------------- data ----------------------------------- #
M = generate_manifold_pointcloud(Plane(); N=20_000)
M = M[:, M[3, :] .> 0]
X = eachcol(M) |> collect

# ----------------------------------- model ---------------------------------- #
K = size(M, 1)
D = 2
inner = 64

encoder = Chain(
    Dense(K=>inner, tanh),
    Dense(inner=>D, tanh),
)

decoder = Chain(
    Dense(D=>inner, tanh),
    Dense(inner=>K, tanh),
)

vae = (decoder ∘ encoder)

# --------------------------------- training --------------------------------- #
"""
Good resource: https://loving-bohr-270bc1.netlify.app/tutorials/2020/11/18/vae_mnist.html
"""


opt = Flux.ADAM(0.0001)

ℓ(x) = (vae(x)   .- x).^2 |> sum
θ = Flux.params(encoder, decoder)

L = []
@track for i in 1:2500
    Flux.Optimise.train!(
        ℓ, θ, X, opt
    )

    l = mean(
        ℓ.(X)
    )
    push!(L, l)
    l < 0.001 && break
end
plot(L, ylim=[0, maximum(L)+.2]) |> display


# -------------------------- plot reconstructed data ------------------------- #
p = plot()
scatter!(eachrow(M)..., ms=5, color=:black, alpha=.5)
Y = hcat(vae.(X)...)
scatter!(eachrow(Y)..., color=:red, 
    xlim=[-2, 2], ylim=[-2, 2], zlim=[-2, 2],
    alpha=.25
)
anim = @animate for x in 0:5:30
    plot!(camera=(x, x))
end
gif(anim, fps=5) |> display




# --------------------------- plot base coordinates -------------------------- #

p2s = plot(title="Encoded coordinates")

# plot encoded coordinates
p1 = plot()
X = range(-1, 1, length=10)
for x ∈ X, y ∈ X, z ∈ X
    encoded = encoder([x, y, z])
    scatter!([encoded[1]], [encoded[2]], marker_z=y, ms=3, label=nothing)
end

# plot "coordinates"
p2 = plot(title="Global coordinates grid")
for v in range(-1, 1, length=15)
    x = range(-1, 1, length=100) |> collect
    y = ones(100) .* v |> collect
    ŷ₁ = hcat(decoder.(eachrow(hcat(x, y)))...)
    ŷ₂  = hcat(decoder.(eachrow(hcat(y,x)))...)
    scatter!(eachrow(ŷ₁)..., alpha=.5, label=nothing, color=:red) 
    scatter!(eachrow(ŷ₂)..., alpha=.5, label=nothing, color=:black)
end
plot(p1, p2, plot(p2, camera=(0, 90)), plot(p2, camera=(90, 0)), size=(800, 800)) |> display





