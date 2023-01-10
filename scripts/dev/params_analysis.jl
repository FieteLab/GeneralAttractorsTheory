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


title = "Vᵢ/Vᵦ - soft RELU"

fld_name = "ring_grid_search"
mfld = "ring"
# B = range(4, 6, length = 2) |> collect
B = [1]
D = range(0.05, 0.5, length = 35) |> collect
V = range(0.05, 0.5, length = 4) |> collect


params = product(B, D, V) |> collect
# colors = getfield.(Palette(indigo, salmon_dark; N=length(B)).colors, :string)
colors = [salmon_dark]

dcolors = getfield.(Palette(indigo, salmon_dark; N = length(D)).colors, :string)
vcolors = getfield.(Palette(indigo, salmon_dark; N = length(V)).colors, :string)

ms = 20

# ----------------------------------- utils ---------------------------------- #
if mfld == "torus"
    can = CAN("torus", cover, n, ξ_t, d_t, k_t; offset_size = 0.1)
elseif mfld == "ring"
    can = CAN("torus", cover, n, ξ_r, d_r, k_r; offset_size = 0.1)
else
    error()
end

# --------------------------------- load data -------------------------------- #
if LOAD
    data = Dict{Symbol,Vector}(
        :δ => [],
        :b => [],
        :v => [],
        :s => [],  # on mfld bump speed
    )


    for δ in D, b in B, v in V
        name = "v_$(v)_δ_$(δ)_b_$(b)_$(mfld)_history"
        s = get_bump_speed(can, fld_name, name)

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
p1 = plot(
    xlabel = "Offset size δ",
    ylabel = "bump v",
    # aspect_ratio = :equal,
    # title = title * " b₀ = $(round(B[1], digits=1))",
)
p2 = plot(
    xlabel = "Offset size δ",
    ylabel = "bump v",
    # aspect_ratio = :equal,
    # title = title * " b₀ = $(round(B[2], digits=1))",
)
# vmin, vmax = minimum(data.v), maximum(data.v)
# plot!(
#     p1,
#     [vmin, vmax],
#     [vmin, vmax],
#     lw = 2,
#     color = "black",
#     alpha = 0.5,
#     ls = :dashdotdot,
#     label = "ideal",
# )
# plot!(
#     p2,
#     [vmin, vmax],
#     [vmin, vmax],
#     lw = 2,
#     color = "black",
#     alpha = 0.5,
#     ls = :dashdotdot,
#     label = nothing,
# )

for (color, v) in zip(vcolors, V)
    for (j, b) in enumerate(B)
        _data = data[(data.v.==v).&(data.b.==b), :]
        p = j == 1 ? p1 : p2
        plot!(
            p,
            _data.δ,
            _data.s,
            lw = 2,
            color = color,
            label = j == 1 ? "v=$(round(v, digits=3))" : nothing,
        )
        plot!(
            [minimum(_data.δ), maximum(_data.δ)],
            [minimum(_data.s), maximum(_data.s)],
            lw = 2,
            color = color,
            ls = :dash,
            alpha = 0.4,
            label = nothing,
        )
    end
end
plot(p1, size = (800, 800)) |> display


p1 = plot(
    xlabel = "Offset size δ",
    ylabel = "max s/ min s",
    # aspect_ratio = :equal,
    # title = title * " b₀ = $(round(B[1], digits=1))",
)
vmin, vmax = V[1], V[2]
for (color, d) in zip(dcolors, D)
    for (j, b) in enumerate(B)
        _data = data[(data.δ.==d).&(data.b.==b), :]

        s_min = _data[_data.v.==vmin, :].s[1]
        s_max = _data[_data.v.==vmax, :].s[1]


        scatter!(p1, [d], [s_max / s_min], lw = 2, color = color, label = nothing)
    end
end


# plot(p1, p2, size = (1000, 800)) |> display
plot(p1, size = (800, 800)) |> display


# # ------------------------ plot v/s for different b\_0 ----------------------- #
# """
# For each δ plot a line showing bump speed over speed for ecah b₀
# """

# plots = []
# for δ in D[1:end-1]
#     plt = plot(
#         title = "δ = $(round(δ; digits=2))",
#         xlabel = "velocity input",
#         ylabel = "bump speed",
#         legend = :topleft,
#         grid = false,
#         ylim = [0, maximum(data.s) * 1.1],
#         # aspect_ratio=:equal,
#     )
#     for (b, color) in zip(B[1:2:end], colors)
#         _data = data[(data.b.==b).&(data.δ.==δ), :]

#         plot!(plt, _data.v, _data.s, lw = 2, color = color, label = "b₀ = $b")
#         scatter!(
#             plt,
#             _data.v,
#             _data.s,
#             ms = 3,
#             color = "white",
#             msc = color,
#             label = nothing,
#         )

#         plot!(
#             _data.v,
#             _data.v,
#             lw = 2,
#             color = :black,
#             alpha = 0.2,
#             ls = :dash,
#             label = nothing,
#         )

#     end
#     push!(plots, plt)
# end
# plot(plots...; size = (1000, 1000)) |> display


# # ----------------------------- plot s/v heatmap ----------------------------- #

# """
# For each v plot a heatmap showing s/v for each δ/b₀
# """

# plots = []
# for v in V
#     plt = plot(
#         title = "v = $(round(v; digits=2))",
#         xlabel = "b₀",
#         ylabel = "δ",
#         grid = false,
#         # ylim=[0, maximum(data.s)*1.1]
#         # aspect_ratio=:equal,
#     )

#     _data = data[data.v.==v, :]
#     scatter!(
#         plt,
#         _data.b,
#         _data.δ,
#         marker_z = _data.s ./ (v * 1),
#         msa = 0,
#         msw = 0,
#         # clims=(0.0, 2.0),
#         clim = (0.02, 0.12),
#         cbar_title = "bump v / input v",
#         color = :bwr,
#         ms = ms,
#         label = nothing,
#     )
#     push!(plots, plt)
# end
# plot(plots...; size = (800, 800)) |> display
