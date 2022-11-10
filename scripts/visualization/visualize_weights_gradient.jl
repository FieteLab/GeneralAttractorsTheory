using Plots
using LinearAlgebra
using ForwardDiff
using Interpolations


using GeneralAttractors.Simulations
using GeneralAttractors.Kernels
using GeneralAttractors.ManifoldUtils
using GeneralAttractors: lerp
import GeneralAttractors.ManifoldUtils: sphere_embedding
using Distances

# --------------------------------- make net --------------------------------- #
n = (64, 64)
function ξ_s(i::Int, j::Int)::Vector
    [lerp(i, n[1], -π, π), lerp(j, n[2], -π / 2, π / 2)]
end
d_s = SphericalAngle()
k_s = DiffOfExpKernel(; λ = 0.75)

cover = CoverSpace(S², S², (x, y) -> [x, y])
can = CAN(
    "sphere",
    cover,
    n,
    ξ_s,
    d_s,
    k_s;
    offset_size = 0.2,
    φ = sphere_embedding,
    offsets = [[1, 1], [1, -1], [-1, 1], [-1, -1]],
)


# ----------------------------------- plot ----------------------------------- #
copyidx = 3


x = [x[1] for x in ξ_s.(1:64, 1)]
y = [x[2] for x in ξ_s.(1, 1:64)]

plt = plot(
    xlim = (minimum(x), maximum(x)),
    ylim = (minimum(y), maximum(y)),
    aspect_ratio = :equal,
)

for (idx, color) in zip([3900, 2500], (:black, :red))
    z = reshape(can.Ws[copyidx][:, idx], (64, 64))

    itp = interpolate((x, y), z, Gridded(Linear()))
    grad = gradient.(Ref(itp), x, y')

    for (i, x) in enumerate(x), (j, y) in enumerate(y)
        Δ = grad[i, j] .* 1.5
        norm(Δ) < 0.01 && continue
        plot!([x, x + Δ[1]], [y, y + Δ[2]], lw = 1, color = color, label = nothing)
        scatter!([x], [y], lw = 1, color = color, label = nothing, msa = 0, msw = 0, ms = 1)
    end
    scatter!([can.X[1, idx]], [can.X[2, idx]], ms = 5, color = color, label = nothing)
end
plt
