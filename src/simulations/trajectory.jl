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
function add_initial_still_phase(X::Matrix, X̄::Matrix, V::Matrix, still, x₀)
    standing = zeros(still, length(x₀))
    for i = 1:length(x₀)
        standing[:, i] .= x₀[i]
    end

    X = vcat(standing, X)
    X̄ = vcat(standing, X̄)
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



"""
    get_closest_neuron(x, X)

Given a point x in a manifold get the coordinates
x̂ ∈ X of the closest neuron given a manifold
metric `d`.
"""
function get_closest_neuron(x, X, d)
    Δ = map(x̂ -> d(x̂, x), eachcol(X))
    return X[:, argmin(Δ)]
end

"""
Ensure a vector `v` has magnitude ∈ [-vmax, vmax]
"""
function enforce_vmax(v, vmax)
    μ = norm(v)
    return if μ > vmax
        return v ./ μ .* vmax
    else
        v
    end
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
    X̄::Matrix   # on manifold trajectory
    V::Matrix
    still::Int # number of warmup frames at the start
end


Trajectory(can::AbstractCAN, args...; kwargs...) = Trajectory(can, can.C.M, args...; kwargs...)



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
        can::AbstractCAN,
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
        smoothing_window=11,
    )   
    
    ψs::Vector = can.C.M.ψs # get manifold vector fields
    n_vfields = length(ψs)

    # get manifold dimensionality
    d = length(M.xmin)

    # get velocity params
    σv = σv isa Number ? repeat([σv], n_vfields) : σv
    μv = μv isa Number ? repeat([μv], n_vfields) : μv
    @assert length(σv) == n_vfields
    @assert length(μv) == n_vfields
    # @assert n_vfields >= d


    # get starting point
    x₀ = x₀ isa Number ? repeat([x₀], d) : x₀
    x₀ = !isnothing(x₀) ? x₀ : rand(M)
    @assert length(x₀) == d "Got x₀: $(x₀) and d=$d"


    # get velocity vector magnitude in components at each frame
    Vs = Vector[]  # magnitude trace along each dimension
    for i = 1:n_vfields
        v = if modality == :piecewise
            piecewise_linear(T, n_piecewise_segments, -vmax:(vmax/100):vmax)
        elseif modality == :constant
            ones(T) .* μv[i]
        else
            x = random_variable(T, μv[i], σv[i]; smoothing_window = smoothing_window) * scale
            clamp!(x, -vmax, vmax)
            ramp = [range(0, 1, length=100)..., ones(T-100)...]
            x .* ramp
        end
        push!(Vs, v)
    end

    """
        Get the sum of the vector fields ψᵢ at position p ∈ M
        and at time step t
    """
    ∑ψ(p, t) = map(i -> Vs[i][t] * ψs[i](p), 1:n_vfields) |> sum


    # first, generate a trajectory on the M mfld
    X, V = Matrix(reshape(zeros(T, d), T, d)), Matrix(reshape(zeros(T, d), T, d))
    X[1, :] =  x₀
    for t in 2:T
        x = X[t-1, :]

        # get a velocity vector
        v = ∑ψ(x, t)

        # make sure vmax magnitude is in range
        # v = enforce_vmax(v, vmax)

        # update on mfld position
        x̂ = x + (v * dt)
        x̂, v_correction = apply_boundary_conditions!(x̂, can.C.M)
        
        # store
        X[t, :] = x̂
        V[t-1, :] = v .* v_correction
    end

    # if the M and N manifolds are the same, we're done
    if can.C.M == can.C.N
        X̄ = X
    else
        # reconstruct the N mfld trajecotry using the cover map ρ
        X̄ = by_column(can.C.ρ, Matrix(X'))' |> Matrix
    end

    # add a still phase
    still > 0 && begin
        X, V = add_initial_still_phase(X, X̄, V, still, x₀)
    end
    return Trajectory(M, X, X̄, V, still)
end


