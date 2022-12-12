import Term.Repr: @with_repr
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
    switches = (1:n_intervals) .* (T / n_intervals) .- (T / n_intervals) |> collect

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
    for i = 1:length(x₀)
        standing[:, i] .= x₀[i]
    end

    X = vcat(standing, X)
    V = vcat(zeros(still, length(x₀)), V)
    return X, V
end


"""
    random_variable(N::Int, μ::Number, σ::Number; smoothing_window=nothing)

1d white noise trace of length N, std `σ` and mean `μ`. Optionally smoothed.
"""
function random_variable(N::Int, μ::Number, σ::Number; smoothing_window = nothing)
    x = (rand(N) .- 0.5) .* σ .+ μ
    isnothing(smoothing_window) || (x = moving_average(x, smoothing_window))
    return x
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
@with_repr struct Trajectory
    M::AbstractManifold
    X::Matrix
    V::Matrix
    still::Int # number of warmup frames at the start
end


Trajectory(can::AbstractCAN, args...; kwargs...) = Trajectory(can.C.M, args...; kwargs...)



"""
Trajectory(
    M::AbstractManifold,
    T::Int = 250,
    dt::Float64 = 0.5,
    σv::Vector = [1, 1, 1],  # variability of each speed field
    μv::Vector = [0, 0, 0],  # bias of each speed field
    x₀ = nothing,
    still = 0,
    vmax = 0.075,
    modality = :piecewise,
    n_piecewise_segments = 3,
)

Generic trajectory constructor. 
Given a list of vector fields (ψᵢ) on a maniofold M,
create a list of random white noise (smoothed) magnitude
trajectories vᵢ for each vec field at each time step.
From these integrate to reconstruct trajectory.
Optionally, the trajectory can be made out of 
segments of constant magnitude traces (modality=:piecewise)
or it can be constant throughout (modality=:constant)

Arguments:
 - M: manifold
 - ψs: vector fields (defined as functions), one for each dimension of the mfld, at least.
 - T, dt: total number of steps and Δt
 - σv: std of the magnitude trace of each respective vfield
 - μv: mean of each magnitude trace
 - still: initial period with no speed
 - x₀: initial condition
 - vmax: cap on magnitude traces of each components of the velocity vector.
 - modality: trajectory modality type
 - scale: used to scale the velocity vectors.

"""
function Trajectory(
    M::AbstractManifold;
    T::Int = 250,
    dt::Float64 = 0.5,
    σv::Union{Number,Vector} = 1,  # variability of each speed field
    μv::Union{Number,Vector} = 0,  # bias of each speed field
    x₀::Union{Nothing,Number,Vector} = nothing,
    still::Int = 0,
    vmax::Number = 0.075,
    modality::Symbol = :random,
    n_piecewise_segments::Int = 3,
    scale::Number = 1,
)
    ψs::Vector = M.ψs # get manifold vector fields

    # get manifold dimensionality
    d = length(M.xmin)

    # get velocity params
    σv = σv isa Number ? repeat([σv], d) : σv
    μv = μv isa Number ? repeat([μv], d) : μv
    @assert length(σv) == d
    @assert length(μv) == d
    @assert length(ψs) >= d


    # get starting point
    x₀ = x₀ isa Number ? repeat([x₀], d) : x₀
    x₀ = !isnothing(x₀) ? x₀ : rand(M)
    @assert length(x₀) == d


    # get velocity vector magnitude in components at each frame
    Vs = Vector[]  # magnitude trace along each dimension
    for i = 1:d
        v = if modality == :piecewise
            piecewise_linear(T, n_piecewise_segments, -vmax:(vmax/100):vmax)
        elseif modality == :constant
            ones(T) .* μv[i]
        else
            random_variable(T, μv[i], σv[i]; smoothing_window = 101) # |> cumsum
        end
        push!(Vs, v)
    end

    """
        Get the sum of the vector fields ψᵢ at position p ∈ M
        and at time step t
    """
    ∑ψ(p, t) = map(i -> Vs[i][t] * ψs[i](p), 1:d) |> sum

    # get position and velocity vectors
    X, V = Matrix(reshape(zeros(T, d), T, d)), Matrix(reshape(zeros(T, d), T, d))
    X[1, :] = x₀
    for t = 1:T
        # get velocity vector `v` at time `t` given vfields at position `x`.
        x = t == 1 ? x₀ : X[t-1, :]
        v = ∑ψ(x, t) .* scale
        v = norm(v) > vmax ? v / norm(v) * vmax : v

        # get position at next time step
        x̂ = x + v * dt
        x̂ = apply_boundary_conditions!(x̂, M)

        # store results
        X[t, :] = x̂
        V[t, :] = v
    end


    # add a still phase
    still > 0 && begin
        X, V = add_initial_still_phase(X, V, still, x₀)
    end
    return Trajectory(M, X, V, still)
end
