using GLMakie
import ForwardDiff: jacobian
import GeneralAttractors: MobiusEuclidean

t = range(-1/2, 1/2, length=20)
θ = range(0, 2π-.1, length=50)

f(t, θ) = [
    (1-t*sin(θ/2))*cos(θ),
    (1-t*sin(θ/2))*sin(θ),
    t*cos(θ/2)
]
f(p) = f(p...)



M = hcat([[x, y] for x in t for y in θ]...)
N = hcat(f.(eachcol(M))...)

dist = MobiusEuclidean(2π)
x = [1/2, 0]
d = map(
    i -> dist(x, M[:, i]), 1:size(M, 2)
)


fig = scatter(eachrow(N)..., color=:black)

# TODO define other vfields

Δx = [0, 1]
Δy(t, θ) = [(-0.1)-t, 0]
Δy2(t, θ) = [0.1-t, 0]
J(t, θ) = jacobian(f, [t, θ])

tovec(x) = [[x] for x in x]
for t in -1/2:.1:1/2, θ in 0:.75:(2π-.75)
    v = J(t, θ) * Δx .* 0.1
    arrows!(tovec(f(t, θ))..., tovec(v)..., color=:green, arrowsize=.05)

    v = J(t, θ) * Δy(t, θ) .* 0.1
    arrows!(tovec(f(t, θ))..., tovec(v)..., color=:red, arrowsize=.05)

    v = J(t, θ-.05) * Δy2(t, θ-.05) .* 0.1
    arrows!(tovec(f(t, θ-.05))..., tovec(v)..., color=:blue, arrowsize=.05)
end

fig





