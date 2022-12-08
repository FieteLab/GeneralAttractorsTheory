using Plots
using Term
install_term_stacktrace()
using DataFrames
using Statistics
import MyterialColors: indigo, salmon_dark, Palette

using GeneralAttractors
using GeneralAttractors.Simulations
import GeneralAttractors.Analysis: get_bump_speed

include("../networks/ring.jl")

SIMULATE = true
fld_name = "calibrate_ring"

dt = 0.5
duration = 200
still = 50  # initialization period        
nframes = (Int ∘ round)(duration / dt)

# θ̇₀ = .15
θ₀ = 3.14 # initialize state at center of mfld
d = ringcan.metric.(θ₀, ringcan.X[1, :])
activate = zeros(length(d))
activate[d.<2.5] .= 1

vmax = 0.2
A = range(0.75, 1, step=.05) |> collect
T = range(0.00, vmax, step=0.05) |> collect # baseline speed
colors = getfield.(Palette(indigo, salmon_dark; N=length(T)).colors, :string)
colors2 = getfield.(Palette(indigo, salmon_dark; N=length(A)).colors, :string)


# --------------------------------- simulate --------------------------------- #


function sim()
    num = 0
    for (i, a) in enumerate(A)
        for θ̇₀ in T
            num += 1
            println(Panel("Running sim $num/$(length(A)*length(T))", style = "red", justify = :center))

            can = CAN(
                "ring",
                cover,
                n,
                ξ_r,
                d_r,
                k_r;
                offset_size = 1.0,
                σ = :softrelu,
                α=a,
            )
            
            trajectory = Trajectory(
                can;
                T = nframes,
                dt = dt,
                σθ = 0.0,
                θ̇₀ = θ̇₀,
                θ₀ = θ₀,
                still = still,
                vmax=vmax
            )
            simulation = Simulation(can, trajectory; η = 0, b₀ = 0.5, τ=5.0)
            
            run_simulation(
                simulation;
                frame_every_n = nothing,
                discard_first_ms = still,
                average_over_ms = 0,
                fps = 10,
                s₀ = 1.0 .* activate,
                savefolder = fld_name,
                savename = "a_$(a)_t_$(θ̇₀)",
        )
        end
    end
end

SIMULATE && sim()



# --------------------------------- analysis --------------------------------- #


# --------------------------------- load data -------------------------------- #
data = Dict{Symbol, Vector}(
    :t=>[],
    :a=>[],
    :s=>[],  # on mfld bump speed
)


for a in A
    for θ̇₀ in T
        name = "a_$(a)_t_$(θ̇₀)_ring_history"
        s = get_bump_speed(ringcan, fld_name, name)

        push!(data[:t], θ̇₀)
        push!(data[:a], a)
        push!(data[:s], s)
    end
end

data = DataFrame(data)
# plt = plot(xlabel="α", ylabel="speed", ylim=[0, .3])
# for (color, θ̇₀) in zip(colors, T)
#     _data = data[data.t .== θ̇₀, :]
#     plot!(_data.a, _data.s, label="θ̇₀: "*string(θ̇₀), lw=2, color=color) 


#     th = findfirst(_data.s .>= θ̇₀)
#     isnothing(th) && continue
#     # hline!([θ̇₀], label=nothing, ls=:dash, color=color)
#     plot!([_data.a[th], _data.a[th]], [0, _data.s[th]], lw=2, ls=:dash, color=color, label=nothing)
#     plot!([_data.a[1], _data.a[th]], [θ̇₀, θ̇₀], lw=2, ls=:dash, color=color, label=nothing)
# end
# display(plt)

plt = plot(xlabel="v", ylabel="speed", ylim=[0, vmax], xlim=[0, vmax], aspect_ratio=:equal,
    legend=:topleft, title="bump speed vs input speed")

vmin, vmax = minimum(data.t), maximum(data.t)
plot!(
    [vmin, vmax], [vmin, vmax], lw=4, color=:black, alpha=.3, ls=:dash, label=nothing
)
for (color, a) in zip(colors2, A)
    _data = data[data.a .== a, :]
    plot!(_data.t, _data.s, label="a: "*string(a), lw=2, color=color) 

end
display(plt)