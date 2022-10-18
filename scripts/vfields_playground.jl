using GLMakie
import Base.Iterators: product


# --------------------------- visualize streamlines -------------------------- #
# GLMakie.inline!(false)

# fig = Figure(resolution=(1200, 1200), viewmode = :fitzoom)
# ax = LScene(fig[1, 1], 
#         scenewk=(; padding=(0, 0,0))
#     )


# # visualize streamplot
# ψ(x, y, z) = Point3f0(z*y, -z*x, -z) 
# ψ̂(x, y, z) = abs(x^2+y^2-1) < 0.2 ? ψ(x, y, z) : Point3f0(0, 0, 0)
# ψ̂(x) = ψ̂(x...)
# streamplot!(
#     ψ̂,
#     -1.2..1.2, -1.2..1.2, -1..1,
#     gridsize = (10, 10, 25), arrow_size = 0.1, linewidth=5
# )

# #visualize circle
# θ = 0:.1:2π
# lines!(sin.(θ), cos.(θ), 0.0 .* θ, linewidth=10, color=:black)
# fig


# --------------------------------- simulate --------------------------------- #
import GeneralAttractors: moving_average
GLMakie.inline!(true)


N = 2500  # number of frames
dt = 0.05

# generate ground trouth data
θ̇ = moving_average(
    rand(N) .- 0.5, 20
) 


θ₀ = 0
θ = θ₀ .+ cumsum(θ̇) .* dt
lines(θ)

# simulate
u =  [0 0 1] .* θ̇  # N × 3
p₀ = [1, 0, 0]

α, β, γ = 5, 5, 5
ṗ(x, y, z) = [-α*z*y, β*z*x, -γ*z]

T = zeros(3, N)
T[:, 1] = p₀ .+ dt .* u[1, :]
for i in 2:N
    δt = ṗ(T[:, i-1]...) + u[i, :]
    T[:, i] = T[:, i-1] + δt*dt


end


f = Figure()

ax = f[1, 1] = Axis(f)
θ̂ = atan.(T[2, :] ./ T[1, :]) .+ θ₀


lines!(
    θ, color=:black; label="real θ"
)
lines!(θ̂, color=:red; label="estimated: θ")
f[1, 2] = Legend(f, ax, "Decoded", framevisible = false)
f
