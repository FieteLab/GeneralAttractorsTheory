using Plots
using Term
install_term_stacktrace()
using DataFrames
using Statistics
import MyterialColors: indigo, salmon_dark, Palette

using GeneralAttractors
using GeneralAttractors.Simulations
import GeneralAttractors.Analysis: get_bump_speed

include("../networks/torus.jl")

SIMULATE = false
fld_name = "calibrate_tours"
can = toruscan

dt = 0.5
duration = 200
still = 50  # initialization period        
nframes = (Int ∘ round)(duration / dt)

# init params
x₀ = [3.15, 3.15] # initialize state at center of mfld
d = map(i -> toruscan.metric(x₀, toruscan.X[:, i]), 1:size(toruscan.X, 2))
activate = zeros(length(d))
activate[d.<2.5] .= 1

# ---------------------------------- params ---------------------------------- #
τ = 5.0
vmax = 0.2


A = range(0.5, 1, step=.5) |> collect
V = range(0.00, vmax, step=0.1) |> collect # baseline speed
colors = getfield.(Palette(indigo, salmon_dark; N=length(V)).colors, :string)
colors2 = getfield.(Palette(indigo, salmon_dark; N=length(A)).colors, :string)




# --------------------------------- simulate --------------------------------- #
function sim()
    num = 0
    for (i, a) in enumerate(A)
        for v0 in V
            num += 1
            println(Panel("Running sim $num/$(length(A)*length(V))", style = "red", justify = :center))

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
                vmax=vmax
            )
            simulation = Simulation(can, trajectory; η = 0, b₀ = 0.5, τ=τ)

            run_simulation(
                simulation;
                frame_every_n = nothing,
                discard_first_ms = still,
                average_over_ms = 0,
                fps = 10,
                s₀ = 1.0 .* activate,
                savefolder = fld_name,
                savename = "a_$(a)_v_$(v0)",
        )
        end
    end
end

SIMULATE && sim()



# --------------------------------- analysis --------------------------------- #
data = Dict{Symbol, Vector}(
    :v=>[],
    :a=>[],
    :s=>[],  # on mfld bump speed
)


for a in A
    for v0 in V
        name = "a_$(a)_v_$(v0)_$(can.name)_history"
        s = get_bump_speed(can, fld_name, name)

        push!(data[:v], v0)
        push!(data[:a], a)
        push!(data[:s], s)
    end
end

data = DataFrame(data)

plt = plot(xlabel="v", ylabel="speed", ylim=[0, vmax], xlim=[0, vmax], aspect_ratio=:equal,
    legend=:topleft, title="bump speed vs input speed")

vmin, vmax = minimum(data.v), maximum(data.v)
plot!(
    [vmin, vmax], [vmin, vmax], lw=4, color=:black, alpha=.3, ls=:dash, label=nothing
)
for (color, a) in zip(colors2, A)
    _data = data[data.a .== a, :]
    plot!(_data.v, _data.s, label="a: "*string(a), lw=2, color=color) 

end
display(plt)





# plt = plot(xlabel="α", ylabel="speed", ylim=[0, .3])
# for (color, θ̇₀) in zip(colors, T)
#     _data = data[data.v .== θ̇₀, :]
#     plot!(_data.a, _data.s, label="θ̇₀: "*string(θ̇₀), lw=2, color=color) 


#     th = findfirst(_data.s .>= θ̇₀)
#     isnothing(th) && continue
#     # hline!([θ̇₀], label=nothing, ls=:dash, color=color)
#     plot!([_data.a[th], _data.a[th]], [0, _data.s[th]], lw=2, ls=:dash, color=color, label=nothing)
#     plot!([_data.a[1], _data.a[th]], [θ̇₀, θ̇₀], lw=2, ls=:dash, color=color, label=nothing)
# end
# display(plt)
