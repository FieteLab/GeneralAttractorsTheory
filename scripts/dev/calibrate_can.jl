using Plots
using Term
install_term_stacktrace()
install_term_logger()

using DataFrames
using Statistics
import MyterialColors: indigo, salmon_dark, Palette

using GeneralAttractors
using GeneralAttractors.Simulations
import GeneralAttractors.Analysis: get_bump_speed

include("../networks/ring.jl")

SIMULATE = true
fld_name = "calibrate_rinf"
can = ringcan

dt = 0.5
duration = 100
still = 50  # initialization period        
nframes = (Int ∘ round)(duration / dt)

# init params
if can.name == "torus"
    x₀ = [3.15, 3.15] # initialize state at center of mfld
    d = map(i -> toruscan.metric(x₀, toruscan.X[:, i]), 1:size(toruscan.X, 2))
    activate = zeros(length(d))
    activate[d.<2.5] .= 1
elseif can.name == "ring"
    x₀ = 3.14 # initialize state at center of mfld
    d = can.metric.(x₀, can.X[1, :])
    activate = zeros(length(d))
    activate[d.<2.5] .= 1

end

# ---------------------------------- params ---------------------------------- #
τ = 5.0
vmax = 0.3
Δv = 0.1

A = range(0.4, 0.5, step = 0.1) |> collect
V = range(0, vmax, step = Δv) |> collect # baseline speed
B = [0.1, 1.0]
colors = getfield.(Palette(indigo, salmon_dark; N = length(V)).colors, :string)
colors2 = getfield.(Palette(indigo, salmon_dark; N = length(A)).colors, :string)
colors3 = getfield.(Palette(indigo, salmon_dark; N = length(B)).colors, :string)




# --------------------------------- simulate --------------------------------- #
function sim()
    num = 0
    for a in A, v0 in V, b0 in B
        num += 1
        println(
            Panel(
                "Running sim $num/$(length(A)*length(V)*length(B))",
                style = "red",
                justify = :center,
            ),
        )

        can.α = a

        trajectory = Trajectory(
            can;
            T = nframes,
            dt = dt,
            σθ = 0.0,
            σv = 0.0,
            μv = v0,
            x₀ = x₀,
            still = still,
            vmax = vmax,
        )
        simulation = Simulation(can, trajectory; η = 0, b₀ = b0, τ = τ)

        run_simulation(
            simulation;
            frame_every_n = nothing,
            discard_first_ms = still,
            average_over_ms = 0,
            fps = 10,
            s₀ = 1.0 .* activate,
            savefolder = fld_name,
            savename = "a_$(a)_v_$(v0)_b_$(b0)",
        )
    end
end

SIMULATE && sim()



# --------------------------------- analysis --------------------------------- #
data = Dict{Symbol,Vector}(
    :v => [],
    :a => [],
    :b => [],
    :s => [],  # on mfld bump speed
)


for a in A, v0 in V, b0 in B
    name = "a_$(a)_v_$(v0)_b_$(b0)_$(can.name)_history"
    s = get_bump_speed(can, fld_name, name)

    push!(data[:v], v0)
    push!(data[:b], b0)
    push!(data[:a], a)
    push!(data[:s], s)
end

data = DataFrame(data)

plt = plot(
    xlabel = "v",
    ylabel = "speed",
    ylim = [0, vmax],
    xlim = [0, vmax],
    aspect_ratio = :equal,
    legend = :topleft,
    title = "bump speed vs input speed",
)

vmin, vmax = minimum(data.v), maximum(data.v)
plot!(
    [vmin, vmax],
    [vmin, vmax],
    lw = 4,
    color = :black,
    alpha = 0.3,
    ls = :dash,
    label = nothing,
)
for (color, a) in zip(colors2, A)
    for (j, (ls, b0)) in enumerate(zip((:solid, :dash, :dashdotdot), B))
        _data = data[(data.a.==a).&(data.b.==b0), :]
        plot!(
            _data.v,
            _data.s,
            label = j == 1 ? "a: " * string(a) : nothing,
            lw = 2,
            color = color,
            ls = ls,
        )
    end
end
display(plt)



# plt = plot(xlabel="α", ylabel="speed", ylim=[0, vmax], xlim=[0, vmax], aspect_ratio=:equal)
# for (color, θ̇₀) in zip(colors, V)
#     for (j, (ls, b0)) in enumerate(zip((:solid, :dash, :dashdotdot), B))
#         _data = data[(data.v .== θ̇₀) .& (data.b .== b0), :]
#         plot!(_data.a, _data.s, 
#                 label= j == 1 ? "θ̇₀: "*string(θ̇₀) : nothing, 
#                 lw=2, color=color, ls=ls) 


#         # th = findfirst(_data.s .>= θ̇₀)
#         # isnothing(th) && continue

#         # plot!([_data.a[th], _data.a[th]], [0, _data.s[th]], lw=2, ls=:dash, color=color, label=nothing)
#         # plot!([_data.a[1], _data.a[th]], [θ̇₀, θ̇₀], lw=2, ls=:dash, color=color, label=nothing)
#     end
# end
# display(plt)
