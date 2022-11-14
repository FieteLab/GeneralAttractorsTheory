

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
function Trajectory(
    M::Manifoldℝ²;
    T::Int = 250,
    σv = 0.5,
    μ = 0.1,
    σθ = 0.5,
    θ₀ = nothing,
    x₀ = nothing,
    y₀ = nothing,
)
    v = rand(T) .* σv .+ μ
    v[v.<0] .= 0

    θ₀ = isnothing(θ₀) ? rand(0:0.2:2π) : θ₀
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

"""
Define random trajectories on the sphere (in ℝ³) by:
i. defining three smooth 1-d random functions
ii. use these to take linear combinations of killing fields of the unit sphere
"""
function Trajectory(
    M::Sphere;
    T::Int = 250,
    σ = [1, 1, 1],
    scale=0.01, 
    x₀=nothing,
    still=0,
    vmax=0.1
)
    dt = 100
    T2 = T*dt

    # get starting point
    x₀ = !isnothing(x₀) ? x₀ : begin
        p = rand(-1:.1:1, 3)
        p ./= norm(p)
    end

    # smooth sum of vector fields
    vx = moving_average((rand(T2).-0.5) .* σ[1], 20dt) |> cumsum
    vy = moving_average((rand(T2).-0.5) .* σ[2], 20dt) |> cumsum
    vz = moving_average((rand(T2).-0.5) .* σ[3], 20dt) |> cumsum


    a = (Int ∘ round)(T2/2)
    vz = ones(T2) .* 0.05
    vy = ones(T2) .* 0.1
    # vz[a:end] .*= 0.0
    vy[1:a] .*= 0.0

    # vy = range(0, .1, length=T2) |> collect
    
    
    clamp!(vx, -vmax, vmax)
    clamp!(vy, -vmax, vmax)
    clamp!(vz, -vmax, vmax)


    ∑ψ(p, i) = scale .* (vx[i]*ψx(p...) + vy[i]*ψy(p...) + vz[i]*ψz(p...))
    X, V = zeros(T2, 3), zeros(T2, 3)
    X[1, :] = x₀
    for t in 2:T2
        x = X[t-1, :]
        v = ∑ψ(x, t)
        x̂ =  X[t-1, :] + v/dt
        x̂ ./= norm(x̂)

        X[t, :] = x̂
        V[t, :] = v
    end

    # add a still phase at the beginning
    X, V = X[1:dt:end, :], V[1:dt:end, :]

    if still > 0
        standing = zeros(still, 3)
        standing[:, 1] .= x₀[1]
        standing[:, 2] .= x₀[2]
        standing[:, 3] .= x₀[3]
        X = vcat(standing, X)
        V = vcat(zeros(still, 3), V)
    end


    return Trajectory(M, X, V)
end
