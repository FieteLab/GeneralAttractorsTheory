# ---------------------------------------------------------------------------- #
#                                  EMBEDDINGS                                  #
# ---------------------------------------------------------------------------- #
"""
Classic sphere embedding in ℝ³
"""
function sphere_embedding end

sphere_embedding(p) = sphere_embedding(p...)
function sphere_embedding(lon, lat)
    ls = atan(tan(lat))
    return [cos(ls) * cos(lon), cos(ls) * sin(lon), sin(ls)]
end




"""
Mobius band embedding in ℝ³

Assumes parameters:
    - t ∈ [-1/2, 1/2]
    - θ ∈ [0, 2π]
"""
function mobius_embedding end

mobius_embedding(t, θ) = [
    (1-t*sin(θ/2))*cos(θ),
    (1-t*sin(θ/2))*sin(θ),
    t*cos(θ/2)
]
mobius_embedding(p) = mobius_embedding(p...)

mobius_embedding(x::Matrix) = hcat(map(mobius_embedding, eachcol(x))...)



function torus_embedding end

function torus_embedding(θ₁, θ₂)
    R, r = 1.0, 0.7
    return [(R + r * cos(θ₁)) * cos(θ₂), (R + r * cos(θ₁)) * sin(θ₂), r * sin(θ₁)]
end
torus_embedding(p) = torus_embedding(p...)