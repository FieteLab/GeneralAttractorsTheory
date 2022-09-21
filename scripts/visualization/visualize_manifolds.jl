import GeneralAttractors.ManifoldUtils: 
        TorusEmbedding,
        SphereEmbedding,
        visualize_manifold,
        MobiusEmbedding,
        CylinderEmbedding

using GLMakie



GLMakie.inline!(false)  # set to false to get interactive 3d renderings


mfld = :torus  # choose which manifold to visualize
colorby, cmap = :d, :inferno
if mfld == :torus
    x = range(-π, π, length=25)
    y = range(-π, π, length=50)
    φ = TorusEmbedding()
elseif mfld == :sphere
    x = range(-π, π, length=25)
    y = range(-π/2, π/2, length=25)
    φ = SphereEmbedding()
    colorby=:z
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

visualize_manifold(x, y, φ; colorby=colorby, cmap=cmap, transparency=false)