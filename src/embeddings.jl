# ---------------------------------------------------------------------------- #
#                                  EMBEDDINGS                                  #
# ---------------------------------------------------------------------------- #


function identity_embedding end

identity_embedding(x, y) = [x, y]
identity_embedding(x, y, z) = [x, y, z]
identity_embedding(p) = identity_embedding(p...)


# ---------------------------------------------------------------------------- #
#                                    SPHERE                                    #
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

sphere_embedding(x, y, z) = [x, y, z]


# ---------------------------------------------------------------------------- #
#                                    MOBIUS                                    #
# ---------------------------------------------------------------------------- #

"""
Mobius band embedding in ℝ³

Assumes parameters:
    - t ∈ [-1/2, 1/2]
    - θ ∈ [0, 2π]
"""
function mobius_embedding end

mobius_embedding(t, θ) =
    [(1 - t * sin(θ / 2)) * cos(θ), (1 - t * sin(θ / 2)) * sin(θ), t * cos(θ / 2)]
mobius_embedding(p) = mobius_embedding(p...)

mobius_embedding(x::Matrix) = hcat(map(mobius_embedding, eachcol(x))...)


# ---------------------------------------------------------------------------- #
#                                 KLEIN BOTTLE                                 #
# ---------------------------------------------------------------------------- #

function klein_embedding end

function klein_embedding(u, v)
    R, r = 0.8, 0.4
    u = u * 2π  # Scale u to [0, 2π]
    v = v * 2π  # Scale v to [0, 2π]
    
    if u < π
        x = (R + r * cos(v)) * cos(u)
        y = (R + r * cos(v)) * sin(u)
        z = r * sin(v)
    else
        x = (R + r * cos(v)) * cos(u)
        y = (R + r * cos(v)) * sin(u)
        z = r * sin(v + π)
    end
    
    return [x, y, z]
end
klein_embedding(p) = klein_embedding(p...)

klein_embedding(x::Matrix) = hcat(map(klein_embedding, eachcol(x))...)




# ---------------------------------------------------------------------------- #
#                                     TORUS                                    #
# ---------------------------------------------------------------------------- #

function torus_embedding end

function torus_embedding(θ₁, θ₂)
    R, r = 0.8, 0.4
    return [(R + r * cos(θ₁)) * cos(θ₂), (R + r * cos(θ₁)) * sin(θ₂), r * sin(θ₁)]
end
torus_embedding(p) = torus_embedding(p...)


# ---------------------------------------------------------------------------- #
#                                   CYLINDER                                   #
# ---------------------------------------------------------------------------- #

function cylinder_embedding end

function cylinder_embedding(θ, z)
    R = 0.8
    return [R * cos(θ), R * sin(θ), z]
end
cylinder_embedding(p) = cylinder_embedding(p...)


# ---------------------------------------------------------------------------- #
#                                     PLANE                                    #
# ---------------------------------------------------------------------------- #

function plane_embedding end

plane_embedding(x, y) = [x, y, 0.0]
plane_embedding(p) = plane_embedding(p...)


# ---------------------------------------------------------------------------- #
#                                     RING                                     #
# ---------------------------------------------------------------------------- #

function ring_embedding end

ring_embedding(Θ::Number) = [cos(θ), sin(θ), 0]
ring_embedding(p::AbstractVector) = ring_embedding(p...)


# ---------------------------------------------------------------------------- #
#                                     LINE                                     #
# ---------------------------------------------------------------------------- #

function line_embedding end

line_embedding(x) = [x, 0, 0]
