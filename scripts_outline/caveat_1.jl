using GeneralAttractors
using GeneralAttractors.ManifoldUtils
using Plots
pyplot()
import GeneralAttractors: sphere_embedding

pts = fibonacci_sphere(350)

p1 = scatter3d(eachrow(pts)..., 
        markersize=5.5, color="#726196", 
        
        label=nothing,
        title = "Fibonacci sampling on S² ⊂ ℝ³"
        )



lon = range(0, 2π; length=25)
lat = range(-π/2, π/2; length=25)

pts2 = [sphere_embedding(o, a) for o in lon for a in lat]
pts2 = hcat(pts2...)
p2 = scatter3d(eachrow(pts2)..., 
        markersize=5.5, color="#726196", 
        
        label=nothing,
        title = "Longitude latitude ebedding"
        
        )


plot(p1, p2, 
        layout=(1,2), size=(800, 400), 
        aspect_ratio=:equal,
        showaxis = false,
        axis=nothing,)