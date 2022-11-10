using GLMakie
using DomainSets
using DomainSets: ×
using ForwardDiff
import ForwardDiff: derivative
using Rotations

GLMakie.inline!(false)

tovec(x) = [[isnan(x) ? 0.0 : x] for x in x]

function sample end
function sample(i::Interval, n=100)::AbstractVector
    range(leftendpoint(i), rightendpoint(i), length=n) |> collect
end

function sample(d::Rectangle; n=100)::Vector{Vector}
    Ω = components(d)
    N = n isa Int ? repeat([n], length(Ω)) : n
    [sample(ω, n) for (ω, n) in zip(Ω, N)]
end


domain = (-π..π)×(-π/2..π/2)
# domain = (0..π)×(0..2π)


pts = sample(domain; n=(20, 10))
M = hcat([[x, y] for x in pts[1] for y in pts[2]]...)  # 2 × n²

"""
Define and visualize vector fields on the sphere embedded in R3
"""


""" embedding """
# φ(θ, ϕ) = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
function φ(lon, lat) 
    ls = atan(tan(lat)) 
    return [ 
            cos(ls) * cos(lon),
            cos(ls) * sin(lon),
            sin(ls),
    ]
end
φ(x) = φ(x...)




∂x = [1, 0, 0]
∂y = [0, 1, 0]
∂z = [0, 0, 1]


# ------------------------------------ viz ----------------------------------- #
S = hcat(map(φ, eachcol(M))...)  # 3 × n

fig = scatter(S[1, :], S[2, :], S[3, :], markersize=25)


for i in 1:size(M, 2)
    p = M[:, i]
    point = φ(p)
    x, y, z = point

    X = 0.2(z*∂y - y*∂z)
    Y = 0.2(z*∂x - x*∂z)
    Z = 0.2(x*∂y - y*∂x)

    arrows!(tovec(point)..., tovec(X)...; arrowsize = Vec3f0(0.1, 0.1, 0.1), color=:black)
    arrows!(tovec(point)..., tovec(-X)...; arrowsize = Vec3f0(0.1, 0.1, 0.1), color=:red)
    # arrows!(tovec(point)..., tovec(Y)...; arrowsize = Vec3f0(0.1, 0.1, 0.1), color=:red)
    # arrows!(tovec(point)..., tovec(Z)...; arrowsize = Vec3f0(0.1, 0.1, 0.1), color=:green)



end
fig