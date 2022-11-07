

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


Trajectory(can::AbstractCAN, args...; kwargs...) = Trajectory(can.M, args...; kwargs...)


"""
    ComputeTrajectory(; T::Int=250, μ=0.1, θ=0.5)

a random walk of duration T over the manifold
ℝ² with average velocity μ and average angular velocity θ.
The trajectory is computed by assuming that an agent moves 
with a certain velocity (randomly drawn from gaussian) and a 
given angular velocity (randomy drawn) reflecting a change in orientation
"""
function Trajectory(M::Manifoldℝ²; T::Int=250, μ=0.1, θ=0.5)

    v = rand(T) .* μ
    θ = cumsum((rand(T) .- 0.5) .* θ)  # orientation

    vx = v .* cos.(θ)
    vy = v .* sin.(θ) 

    x = cumsum(vx)
    y = cumsum(vy)

    return Trajectory(M, hcat(x, y), hcat(vx, vy))
end





