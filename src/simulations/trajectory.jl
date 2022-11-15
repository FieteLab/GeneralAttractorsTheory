# ----------------------------------- utils ---------------------------------- #
"""
    function piecewise_linear(
        T::Int,  # total number of frames
        n_intervals::Int,  # number of different liear segments
        rng::AbstractVector,  # values ranges
    )

Construct a piece-wise linear 1d trajectory
with `T` many frames, and `n_intervals` values
drawn from an iterable `rgn`.
"""
function piecewise_linear(
    T::Int,  # total number of frames
    n_intervals::Int,  # number of different liear segments
    rng::AbstractVector,  # values ranges
)
    # get values and switches timepoints
    vals = rand(rng, n_intervals)
    switches = (1:n_intervals) .* (T/n_intervals ) .- (T/n_intervals ) |> collect

    function x(t)
        idx = findlast(switches .<= t)
        idx = isnothing(idx) ? 1 : idx
        return vals[idx]
    end


    # get value at each step
    t = 1:T |> collect
    return x.(t)
end


"""
Prepend an initial stationary pahse to a trajectory. 
Used to let the network settle in into a selected
state before starting a simulation with velocity inputs.
"""
function add_initial_still_phase(X::Matrix, V::Matrix, still, x₀)
    standing = zeros(still, length(x₀))
    for i in 1:length(x₀)
        standing[:, i] .= x₀[i]
    end

    X = vcat(standing, X)
    V = vcat(zeros(still, length(x₀)), V)
    return X, V
end



# ---------------------------------------------------------------------------- #
#                                  TRAJECTORY                                  #
# ---------------------------------------------------------------------------- #
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
    still::Int # number of warmup frames at the start
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
    dt::Float64=0.5,
    σv = 0.25,
    μv = 0.05,
    vmax = 0.1,
    σθ = 0.5,
    θ₀ = nothing,
    x₀ = nothing,
    y₀ = nothing,
    still = 100,
)   
    # get speed and orientation
    v = rand(T) .* σv .+ μv
    v[v.<vmax] .= vmax
    v[v.<0] .= 0

    θ₀ = isnothing(θ₀) ? rand(0:0.2:2π) : θ₀
    θ̇ = moving_average(rand(T), 11) .- 0.5
    θ̇ = cumsum(θ̇ .* σθ) .+ θ₀ # orientation

    # get velocity at each component
    vx = v .* cos.(θ̇)
    vy = v .* sin.(θ̇)

    # get random initial position
    x₀ = isnothing(x₀) ? rand(-20:2:20) : x₀
    y₀ = isnothing(y₀) ? rand(-20:2:20) : y₀
    
    # Get trajectory
    X, V = zeros(T, 2), zeros(T, 2)
    X[1, :] = [x₀, y₀]
    V[:, 1] = vx
    V[:, 2] = vy
    for t in 2:T
        X[t, :] = X[t-1, :] + V[t, :]*dt
    end

    # finalize
    still > 0 && begin
        X, V = add_initial_still_phase(X, V, still, X[1, :])
    end
    return Trajectory(M, X, V, still)
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
    x₀=nothing,
    still=0,
    vmax=0.075,
    modality=:piecewise,
    n_piecewise_segments=3,
)
    dt = 100
    T2 = T*dt

    # get starting point
    x₀ = !isnothing(x₀) ? x₀ : begin
        p = rand(-1:.1:1, 3)
        p ./= norm(p)
    end

    # get vfield "activation" at each frame
    if modality == :piecewise
        vx = piecewise_linear(T2, n_piecewise_segments, -vmax:(vmax/100):vmax)
        vy = piecewise_linear(T2, n_piecewise_segments, -vmax:(vmax/100):vmax)
        vz = piecewise_linear(T2, n_piecewise_segments, -vmax:(vmax/100):vmax)
    elseif modality == :constant
        vx = ones(T2) .* σ[1]
        vy = ones(T2) .* σ[2]
        vz = ones(T2) .* σ[3]
    else
        vx = moving_average((rand(T2).-0.5) .* σ[1], 20dt) |> cumsum
        vy = moving_average((rand(T2).-0.5) .* σ[2], 20dt) |> cumsum
        vz = moving_average((rand(T2).-0.5) .* σ[3], 20dt) |> cumsum
    end
    
    clamp!(vx, -vmax, vmax)
    clamp!(vy, -vmax, vmax)
    clamp!(vz, -vmax, vmax)

    # get position and velocity vectors
    ∑ψ(p, i) = vx[i]*ψx(p...) + vy[i]*ψy(p...) + vz[i]*ψz(p...)
    X, V = zeros(T2, 3), zeros(T2, 3)
    X[1, :] = x₀
    for t in 2:T2
        x = X[t-1, :]
        v = ∑ψ(x, t)

        #! REMOVE
        v = v ./ norm(v) .* vmax
        # if t < 100*dt
        #     v .*= t/(100*dt)
        # end

        x̂ =  X[t-1, :] + v/dt
        x̂ ./= norm(x̂)

        X[t, :] = x̂
        V[t, :] = v
    end

    # add a still phase at the beginning
    X, V = X[1:dt:end, :], V[1:dt:end, :]

    still > 0 && begin
        X, V = add_initial_still_phase(X, V, still, x₀)
    end
    return Trajectory(M, X, V, still)
end

