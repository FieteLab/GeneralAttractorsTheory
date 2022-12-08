using Plots

using GeneralAttractors
import GeneralAttractors.Analysis: get_bump_speed
using DataFrames
import MyterialColors: indigo, salmon_dark, Palette
using Statistics

"""
run_sims_grid.jl runs lots of simulations with varying 
parameters. Here we plot stuff
"""

LOAD = true


fld_name = "params_grid_search_softrelu_v"
# B = range(0.01, 2, length=1) |> collect
B = range(4, 6, length=2) |> collect
D = range(0.8, 1.2, length=4) |> collect
V = range(0.2, 0.5, length=20) |> collect


params = product(B, D, V) |> collect
# colors = getfield.(Palette(indigo, salmon_dark; N=length(B)).colors, :string)
colors = [salmon_dark]

dcolors = getfield.(Palette(indigo, salmon_dark; N=length(D)).colors, :string)

ms = 20

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

# ----------------------- plot v/s for different delta ----------------------- #
plt = plot(xlabel="input v", ylabel="bump v", aspect_ratio=:equal)
vmin, vmax = minimum(data.v), maximum(data.v)
plot!([vmin, vmax], [vmin, vmax], lw=2,color="black", alpha=.5, ls=:dashdotdot, label="ideal")
for (color, δ) in zip(dcolors, D)
    for (j, b) in enumerate(B)
        _data = data[(data.δ .== δ) .& (data.b .== b), :]
        plot!(_data.v, _data.s, lw=2, color=color, 
            label=j == 1 ? "δ=$(round(δ, digits=3))" : nothing, ls= j==1 ? :dash : :solid)
    end
end
display(plt)


# ------------------------ plot v/s for different b\_0 ----------------------- #
"""
For each δ plot a line showing bump speed over speed for ecah b₀
"""

# plots = []
# for δ in D[1:end-1]
#     plt = plot(
#             title="δ = $(round(δ; digits=2))",
#             xlabel="velocity input", ylabel="bump speed", 
#             legend=:topleft, 
#             grid=false, 
#             ylim=[0, maximum(data.s)*1.1]
#             # aspect_ratio=:equal,
#             )
#     for (b, color) in zip(B[1:2:end], colors)
#         _data = data[(data.b .== b) .& (data.δ .== δ), :]

#         plot!(plt, _data.v, _data.s, lw=2, color=color, label="b₀ = $b")
#         scatter!(plt, _data.v, _data.s, ms=3, color="white", msc=color, label=nothing)

#         plot!(_data.v, _data.v, lw=2, color=:black, alpha=.2, ls=:dash, label=nothing)

#     end
#     push!(plots, plt)
# end
# plot(plots...; size=(1000, 1000)) |> display


# ----------------------------- plot s/v heatmap ----------------------------- #

"""
For each v plot a heatmap showing s/v for each δ/b₀
"""

# plots = []
# for v in V
#     plt = plot(
#         title="v = $(round(v; digits=2))",
#         xlabel="b₀", ylabel="δ", 
#         grid=false, 
#         # ylim=[0, maximum(data.s)*1.1]
#         # aspect_ratio=:equal,
#         )

#     _data = data[data.v .== v, :]
#     scatter!(plt, 
#             _data.b, _data.δ, 
#             marker_z=_data.s./(v * 1),
#             msa=0, msw=0,
#             # clims=(0.0, 2.0),
#             clim=(0.02, 0.12),
#             cbar_title="bump v / input v",
#             color=:bwr,
#             ms=ms,  label=nothing,
#         )
#     push!(plots, plt)
# end
# plot(plots...; size=(800, 800)) |> display
