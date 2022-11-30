using Plots

using GeneralAttractors
import GeneralAttractors.Analysis.ManifoldAnalysis: population_average
import GeneralAttractors.Simulations: decode_peak_location
using DataFrames
import MyterialColors: indigo, salmon_dark, Palette


"""
run_sims_grid.jl runs lots of simulations with varying 
parameters. Here we plot stuff
"""

LOAD = false

B = range(0.01, 0.5, length=5)
D = range(0.1, 0.25, length=3)
V = range(0.01, 0.3, length=30)
params = product(B, D, V) |> collect
@info "Setting up" length(params)
fld_name = "params_sims_lin"

colors = getfield.(Palette(indigo, salmon_dark; N=length(B)).colors, :string)


# ----------------------------------- utils ---------------------------------- #
can = CAN(
    "torus",
    cover,
    n,
    ξ_t,
    d_t,
    k_t;
    offset_size = 0.1,
    )



function ∑(x)
    n, _, m = size(x)
    real.(reshape(mean(x, dims = 2), (n, m)))
end


function get_bump_speed(sim_name::String)::Float64
    # load state history
    history = load_simulation_history(fld_name, sim_name)
    s = ∑(history.S)[:, 30:end]

    # get peak location speed
    peak_location = hcat(
        map(
            st -> decode_peak_location(st, can), 
            eachcol(s)
            )...
        )

    on_mfld_speed = map(
        i -> can.metric(
            peak_location[:, i], peak_location[:, i-1]
            ), 
        2:size(peak_location, 2)
    )
    average_speed = sum(on_mfld_speed)/length(on_mfld_speed)  # tot displacement over time
    return average_speed
end

# --------------------------------- load data -------------------------------- #
if LOAD
    data = Dict{Symbol, Vector}(
        :δ=>[],
        :b=>[],
        :v=>[],
        :s=>[],  # on mfld bump speed
    )


    for δ in D, b in B, v in V
        name = "v_$(v)_δ_$(δ)_b_$(b)_torus_history"
        s = get_bump_speed(name)

        push!(data[:δ], δ)
        push!(data[:b], b)
        push!(data[:v], v)
        push!(data[:s], s)
    end

    data = DataFrame(data)
    println(data)
end


# ---------------------------------------------------------------------------- #
#                                     plot                                     #
# ---------------------------------------------------------------------------- #

plt = plot(xlabel="velocity input", ylabel="bump speed", 
        legend=:topleft, 
        grid=false, 
        # ylim=[0, 0.1], xlim=[0, 0.1],
        # aspect_ratio=:equal,
        )
for (b, color) in zip(B, colors)
    _data = data[(data.b .== b) .& (data.δ .== D[3]), :]

    plot!(plt, _data.v, _data.s, lw=2, color=color, label="b₀ = $b")
    scatter!(plt, _data.v, _data.s, ms=3, color="white", msc=color, label=nothing)
end

plt
