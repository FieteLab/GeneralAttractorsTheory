using Plots

using GeneralAttractors
import GeneralAttractors.Analysis.ManifoldAnalysis: population_average
import GeneralAttractors.Simulations: decode_peak_location
using DataFrames
import MyterialColors: indigo, salmon_dark, Palette
using Statistics

"""
run_sims_grid.jl runs lots of simulations with varying 
parameters. Here we plot stuff
"""

LOAD = false


fld_name = "params_grid_search_relu"
B = range(.1, 10, length=10) |> collect
D = range(.1, 1.7, length=20) |> collect
V = range(0.05, 0.5, length=4) |> collect
colors = getfield.(Palette(indigo, salmon_dark; N=length(B)).colors, :string)
ms = 12

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
    average_speed = sum(on_mfld_speed)/(length(on_mfld_speed)*dt)  # tot displacement over time
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
end


# ---------------------------------------------------------------------------- #
#                                     plot                                     #
# ---------------------------------------------------------------------------- #

# ------------------------ plot v/s for different b\_0 ----------------------- #
"""
For each δ plot a line showing bump speed over speed for ecah b₀
"""

plots = []
for δ in D[1:end-1]
    plt = plot(
            title="δ = $(round(δ; digits=2))",
            xlabel="velocity input", ylabel="bump speed", 
            legend=:topleft, 
            grid=false, 
            ylim=[0, maximum(data.s)*1.1]
            # aspect_ratio=:equal,
            )
    for (b, color) in zip(B[1:2:end], colors)
        _data = data[(data.b .== b) .& (data.δ .== δ), :]

        plot!(plt, _data.v, _data.s, lw=2, color=color, label="b₀ = $b")
        scatter!(plt, _data.v, _data.s, ms=3, color="white", msc=color, label=nothing)

        plot!(_data.v, _data.v, lw=2, color=:black, alpha=.2, ls=:dash, label=nothing)

    end
    push!(plots, plt)
end
plot(plots...; size=(1000, 1000)) |> display


# ----------------------------- plot s/v heatmap ----------------------------- #

"""
For each v plot a heatmap showing s/v for each δ/b₀
"""

plots = []
for v in V
    plt = plot(
        title="v = $(round(v; digits=2))",
        xlabel="b₀", ylabel="δ", 
        grid=false, 
        # ylim=[0, maximum(data.s)*1.1]
        # aspect_ratio=:equal,
        )

    _data = data[data.v .== v, :]
    scatter!(plt, 
            _data.b, _data.δ, 
            marker_z=_data.s./v,
            msa=0, msw=0,
            clims=(0.0, 2.0),
            color=:bwr,
            ms=ms,  label=nothing,
        )
    push!(plots, plt)
end
plot(plots...; size=(800, 800)) |> display
