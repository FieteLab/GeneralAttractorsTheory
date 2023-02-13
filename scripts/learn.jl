using Flux
using Flux: @epochs, onehotbatch, mse, throttle
using Plots

using GeneralAttractors
using GeneralAttractors.Simulations
import GeneralAttractors.Simulations: decode_peak_location, Decoder


warmup_duration = 1000  
trial_duration = 20

x₀ = [3.14, 0]
b₀ = 1.0
τ = 5.0

# ------------------------------ network params ------------------------------ #
d , can_size = 2, 28^2
n_hidden = 256
n_training_trajectories = 1000

# ---------------------------------- get CAN --------------------------------- #
can = torus_maker(:single; n=28)

# run CAN with no inputs to initialize state at x₀
warmup_trajectory = ConstantTrajectory(
    can;
    T=warmup_duration,
    x₀=x₀,
)

warmup, _ = run_simulation(
    Simulation(can, warmup_trajectory; η=0.0, b₀=b₀);
    discard_first_ms=0,
    average_over_ms=0,
)


# ----------------------------------- train ---------------------------------- #

function make_network()
    encoder = Chain(
        Dense(2d => n_hidden, relu),  # ẋ, x̂ → ṡ
        Dense(n_hidden => n_hidden, relu),
        Dense(n_hidden => can_size, relu),
    )

    decoder = Chain(
        Dense(can_size => n_hidden, relu),
        Dense(n_hidden => n_hidden, relu),
        Dense(n_hidden => 2d, relu),
    )

    model = Chain(encoder, decoder)

    opt_state = Flux.setup(Flux.Adam(0.01), model)
    return model, encoder, decoder, opt_state
end

function loss(network_input, network_output) # , x_t, x̂_t)::Float64
    autoencoder_loss = mse(network_output, network_input)
    # can_loss = mse(x_t, x̂_t)
    return autoencoder_loss # + can_loss
end



function initialize_can()
    # initialize can state at end of warmup
    S = warmup_history.S[:, 1, end]

    # initialize S -> X decoder
    x̂ = decode_peak_location(S, can)
    decoder = Decoder(
        trajectory.X[1, :],
        x̂;
        decoding_offset = trajectory.X[1, :] .- x̂,
    )

    return S, decoder
end


function training_loop()
    # make network
    model, encoder, decoder, opt_state = make_network()

    @info "Got network" model encoder decoder opt_state

    # get trajectory
    # TODO: make this within the loop to change trajectory at each epoch
    trajectory = Trajectory(
        can;
        T=max(trial_duration, 500),  #  needs to be long enough to create a trajectory
        vmax=0.05,
        still=0, 
        x₀=x₀,
    )
    

    # ruin over random trajectories
    for training_run in 1:n_training_trajectories
        # initialize CAN
        S, can_decoder = initialize_can()

        # run training loop over the entire trajectory
        run_loss = 0
        for t in 1:trial_duration
            # prepare inputs to network
            x_t, v_t = trajectory.X[t, :], trajectory.V[t, :]
            network_input = vcat(x_t, v_t)

  
            # get loss
            # run_loss += loss(network_input, network_input, x_t, x̂_t)

            # update network
            grads = Flux.gradient(model) do m
                # encode/decode
                # ω = encoder(network_input)  # input to the CAN
                network_output = m(network_input)  # reconstructed input

                # # update CAN's state given ω
                # ṡ = (can.W * S .+ ω .+ b₀) .|> can.σ
                # S = (ṡ - S)/τ

                # # decode x from can's state
                # x̂_t, _ = can_decoder(S, can)

                l = loss(network_input, network_output) # , x_t, x̂_t)
                run_loss += l
                l
            end

            isnothing(grads[1]) && error("grads are nothing, something's wrong in the training loop.")

            Flux.update!(opt_state, model, grads[1])
        end
        @info "RUN: {bold green}$(training_run){/bold green} loss: {red}$(round(run_loss; digits=3)){/red}"
    end
end

encoder = training_loop()