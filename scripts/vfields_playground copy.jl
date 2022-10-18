using GLMakie
import Base.Iterators: product


# --------------------------- visualize streamlines -------------------------- #
GLMakie.inline!(true)

fig = Figure(resolution=(1200, 1200), viewmode = :fitzoom)
ax = LScene(fig[1, 1], 
        scenewk=(; padding=(0, 0,0))
    )


relu(x) = max(x, 0)
id(x) = x

# # visualize streamplot
r, θ = 0, 1
W = [r θ; -θ r]
b = [1, 1] .* 0
println("Eigen values")
println.(eigen(W).values)


ψ(x::AbstractVector) = id.(W * x + b) |> Point2f0 
# ψ(x::AbstractVector) = W * relu.(x) |> Point2f0 
streamplot!(
    ψ,
    -1.5..1.5, -1.5..1.5,
    gridsize = (40, 40), arrow_size = 0.05, linewidth=5,
    aspect_ratio=:equal
)

fig
