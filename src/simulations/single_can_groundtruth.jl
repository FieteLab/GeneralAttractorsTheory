


"""
initialize CAN state and decoder at end of warmup
"""
function initialize_can_with_warmup(can::SingleCAN, warmup::History, trajectory::Trajectory)
    # initialize can state at end of warmup
    S = warmup.S[:, 1, end]

    # initialize S -> X decoder
    x̂ = decode_peak_location(S, can)
    decoder = Decoder(
        trajectory.X[1, :],
        x̂;
        decoding_offset = [0], # trajectory.X[1, :] .- x̂,
    )

    return S, decoder
end



"""
jacobian of the cover map of a can
"""
function cover_map_jacobian(x::Vector, can::SingleCAN)
    J = jacobian(can.C.ρ, x)
    J[isnan.(J)] .= 0
    return J
end


"""
    get_W_partial_derivatives!

Get the partial derivatives of the weight matrix 
along each direction at neuron `i`.

Arguments:
- `o`: a matrix to store the partial derivatives
- `i`: the index of the neuron to get the partial derivatives at
- `can`: the can to get the partial derivatives from
- `can_n`: the number of neurons in the can
"""
function get_W_partial_derivatives!(o::Matrix, i::Int, can::SingleCAN, can_n::Int)
    o .*= 0
    Δx = can.X[1, :] .- can.X[1, i]
    wx = can.kernel(Δx)
    ∂w∂x = begin
        o[2:end, :] = diff(reshape(wx, (can_n, can_n)); dims=1)
        vcat(o...)
    end

    o .*= 0
    Δy = can.X[2, :] .- can.X[2, i]
    wy = can.kernel(Δy)
    ∂w∂y = begin
        o[:, 2:end] = diff(reshape(wy, (can_n, can_n)); dims=2)
        vcat(o...)
    end

    return ∂w∂x, ∂w∂y
end

function get_W_partial_derivatives!(o::Vector, i::Int, can::SingleCAN)
    o .*= 0
    Δx = can.X[1, :] .- can.X[1, i]
    wx = can.kernel(Δx)
    ∂w∂x = begin
        o[2:end] = diff(wx)
        o
    end
    return ∂w∂x
end



"""
    generate_groundtruth_data(trajectory)

Given a trajectory, generate ground-truth data on what `ω`,
the velocity-dependent input, to the network should be using a multiplictive 
weight scaling the partial derivatives of the weights at the "correct" location. 
The CAN isinitialized at the end of a warmup constant trajectory

to test:
    x̄, _ =  test()


    plot(eachcol(trajectory.X)..., lw=4, color=:black)
    plot!(eachrow(hcat(x̄...))..., lw=2, color=:red)

"""
function generate_groundtruth_data(
        can::SingleCAN, trajectory::Trajectory, warmup::History;
        α::Number = -110, b₀=1, τ=5
    )::Tuple{Vector, Vector, Matrix}
    can_n = can.n[1]
    S, can_decoder = initialize_can_with_warmup(can, warmup, trajectory)    
    o = zeros(can.n)  # to store derivatives of weights

    n_steps = size(trajectory.X, 1)
    history = zeros(length(S), n_steps)

    x̄, Ω = [], []
    for t in 1:n_steps
        i = argmax(S)

        x = trajectory.X[t, :]
        ẋ = trajectory.V[t, :]
        J = cover_map_jacobian(x, can)
        v = J*ẋ

        # get and sum partial derivatives to move the state
        if can.d > 1
            ∂w∂x, ∂w∂y = get_W_partial_derivatives!(o, i, can, can_n)            
            ω = α * v[1] * ∂w∂x + α * v[2] * ∂w∂y
        else
            ∂w∂x = get_W_partial_derivatives!(o, i, can)
            ω = α * v[1] * ∂w∂x
        end

        Ṡ = (can.W * S .+ ω .+ b₀) .|> can.σ
        S += (Ṡ - S)/τ

        push!(x̄, can_decoder(S, can)[1])
        push!(Ω, ω)
        history[:, t] .= S
    end

    return x̄, Ω, history
end
