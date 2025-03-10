

"""
Visualize distances and connectivity on the lattices
of each CAN topology
"""

include("settings.jl")

using Distances
import GeneralAttractors: plot_distance_function, SphericalDistance, MobiusEuclidean


# ----------------------------- plot connectivity ---------------------------- #

# for network in networks
#     if network != "klein"
#         continue
#     end
#     can = make_single_can(network; cy_extent = 1)

#     fig = if can.d == 2
#         W = reshape(can.W[1,:], can.n)' |> Matrix
#         w_x = range(can.X[1,1], can.X[1,end]; length=can.n[1])
#         w_y = range(can.X[2,1], can.X[2,end]; length=can.n[2])
    
#         main_plot = contourf(
#             w_x, w_y,
#             W, title = "Connectivity matrix ($network)",
#             aspect_ratio = :equal,
#             color=:bwr,
#             linewidth = 0.25,
#             msc=:black,
#             lc=:black,
#             grid = false,
#             )

#         top_plot = plot(w_x, W[1, :], label=nothing; 
#                 lw=3, color=:black, grid=false
#             )

            
#         side_plot = plot(
#                 W[:, 1], w_y, label=nothing; 
#                 lw=3, color=:black, grid=false
#         )
        

#         plot(
#             plot(; plot_remove_axes...), 
#             top_plot, 
#             side_plot, 
#             main_plot, 
#             # grid(2, 2, heights=[.2, .8], widts=[.2, .8]),
#             size=(1000, 1000)
#         )
        
#     else
#         n = round(Int, can.n[1]/2)
#         w = can.W[n, :]
#         plot(
#             can.X[1, :], w, 
#             label=nothing, xlabel="θ'", ylabel="W(θ, θ')", title="Connectivity matrix ($network)",
#             lw=2, color=:black, 
#             fillrange = maximum(w),
#             fillalpha = 0.25,
#             fillcolor = :black
#         )
#     end

#     save_plot(supervisor, fig, "f3_A_connectivity_$network")
#     display(fig)
# end



# ------------------------------- plot distance ------------------------------ #


d = Euclidean()
p = plot_distance_function(d, plot_title = "Euclidean distance")
save_plot(supervisor, p, "f3_A__distance_line")

d = PeriodicEuclidean([2π])
p = plot_distance_function(d, plot_title = "Ring distance")
save_plot(supervisor, p, "f3_A__distance_ring")

d = PeriodicEuclidean([2π, 2π])  # distance function over a torus manifold
p = plot_distance_function(d, plot_title = "Torus distance"; points = [[1, 1]])
save_plot(supervisor, p, "f3_A__distance_torus")


d = PeriodicEuclidean([Inf, Inf])  
p = plot_distance_function(d, plot_title = "Plane distance"; points = [[.2, .2]])
save_plot(supervisor, p, "f3_A__distance_plane")

d = PeriodicEuclidean([2π, Inf])
p = plot_distance_function(d, plot_title = "Cylinder distance", layout = (4, 1))
save_plot(supervisor, p, "f3_A__distance_cylinder")

d = MobiusEuclidean()
p = plot_distance_function(
    d,
    plot_title = "Mobius distance",
    layout = (1, 4),
    colorbar = nothing,
)
save_plot(supervisor, p, "f3_A__distance_mobius")


d = KleinBottleEuclidean()
p = plot_distance_function(
    d,
    plot_title = "Klein Bottle distance",
    layout = (2, 2),
    colorbar = nothing,
)
save_plot(supervisor, p, "f3_A__distance_klein")

d = SphericalDistance()
p = plot_distance_function(
    d,
    plot_title = "Spherical (angle) distance",
    layout = (2, 2),
    size = (800, 600),
)
save_plot(supervisor, p, "f3_A__distance_spherical")
