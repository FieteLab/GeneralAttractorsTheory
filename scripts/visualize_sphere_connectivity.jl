using GeneralAttractors
using GLMakie


include("./networks/sphere.jl")


idx = 2000
W = spherecan.Ws[1][idx, :]

plt = GLMakie.scatter(eachrow(spherecan.X)..., markersize=50, color=W, alpha=0.1)

for i in 1:100:size(spherecan.X, 2)
    x = spherecan.X[:, i]
    
    for (j, color) in zip((1, 3, 6),(:red, :green, :blue))
        v = spherecan.offsets[j].ψ(x) .* 0.1
        arrows!([[y] for y in x]..., [[y] for y in v]..., color=color)
    end
end

println("done")
plt