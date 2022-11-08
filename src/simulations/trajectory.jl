

"""
struct Trajectory

Stores information of a trajectory on a manifold `M`.
The trajectory is defined by a matrix V (T×d) of velocity
along each of d-many dimensions at each time step in T. 
X (T×d) is the coordinates of the trajectory's trace.
"""
struct Trajectory
    M::AbstractManifold
    X::Matrix
    V::Matrix
end


Trajectory(can::AbstractCAN, args...; kwargs...) = Trajectory(can.C.M, args...; kwargs...)


"""
    ComputeTrajectory(; T::Int=250, μ=0.1, θ=0.5)

a random walk of duration T-many steps over the manifold
ℝ² with average velocity μ and average angular velocity θ.
The trajectory is computed by assuming that an agent moves 
with a certain velocity (randomly drawn from gaussian) and a 
given angular velocity (randomy drawn) reflecting a change in orientation
"""
function Trajectory(M::Manifoldℝ²; T::Int = 250, σv = 0.5, μ = 0.1, σθ = 0.5, θ₀=nothing, x₀=nothing, y₀=nothing)
    v = rand(T) .* σv .+ μ

    θ₀ = isnothing(θ₀) ? rand(0:.2:2π) : θ₀
    θ̇ = moving_average(rand(T), 11) .- 0.5
    θ̇ = cumsum(θ̇ .* σθ) .+ θ₀ # orientation

    vx = v .* cos.(θ̇)
    vy = v .* sin.(θ̇)

    x₀ = isnothing(x₀) ? rand(-20:2:20) : x₀
    y₀ = isnothing(y₀) ? rand(-20:2:20) : y₀
    x = cumsum(vx) .+ x₀
    y = cumsum(vy) .+ y₀

    return Trajectory(M, hcat(x, y), hcat(vx, vy))
end
