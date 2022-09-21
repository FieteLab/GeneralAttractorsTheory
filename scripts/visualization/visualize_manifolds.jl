import GeneralAttractors.ManifoldUtils: 
        TorusEmbedding,
        SphereEmbedding,
        visualize_manifold,
        MobiusEmbedding,
        CylinderEmbedding

using GLMakie


GLMakie.activate!()
GLMakie.inline!(true)  # set to false to get interactive 3d renderings


mfld = :cylinder  # choose which manifold to visualize


colorby, cmap = nothing, nothing
if mfld == :torus
    x = range(-π, π, length=100)
    y = range(-π, π, length=100)
    φ = TorusEmbedding()
elseif mfld == :sphere
    x = range(-π, π, length=100)
    y = range(-π/2, π/2, length=100)
    φ = SphereEmbedding()
elseif mfld == :mobius
    x = range(-π, π, length=100)
    y = range(-1, 1, length=5)
    φ = MobiusEmbedding()
    colorby, cmap = :z, :inferno
elseif mfld == :cylinder
    x = range(-π, π, length=100)
    y = range(-1, 1, length=5)
    φ = CylinderEmbedding()
    colorby, cmap = :z, :inferno
else
    error("Manifold name not recognized: $mfld")
end

visualize_manifold(x, y, φ; colorby=colorby, cmap=cmap)