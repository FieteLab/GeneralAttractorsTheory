"""
Collection of code useful to visualize manifolds
"""
module ManifoldUtils
    using GLMakie
    using Parameters

    using Colors
    


    GLMakie.activate!()

    export TorusEmbedding, SphereEmbedding, visualize_manifold, MobiusEmbedding, CylinderEmbedding

    # ---------------------------------------------------------------------------- #
    #                                  EMBEDDINGS                                  #
    # ---------------------------------------------------------------------------- #
    """
        Standard embedding functions for 1-2d manifolds in ℝ³
    """
    abstract type AbstractEmbedding end
    """ embed a single point """
    (e::AbstractEmbedding)(p::Vector) = e.φ(p)

    """ get embedding coordinates of all points in a manifold """
    function (e::AbstractEmbedding)(x::AbstractVector, y::AbstractVector)::Tuple{Matrix, Matrix, Matrix}
        pts = [
            e.φ([x, y]) for x in x, y in y
        ]
        X = [p[1] for p in pts]
        Y = [p[2] for p in pts]
        Z = [p[3] for p in pts]
        return X, Y, Z
    end


    @with_kw struct CylinderEmbedding <: AbstractEmbedding
        φ::Function = (p) -> begin
            x, y = p
            return [
                cos(x),
                sin(x),
                y
            ]
        end
    end



    @with_kw struct TorusEmbedding <: AbstractEmbedding
        φ::Function = (p) -> begin
            R, r = 1.0, 0.5
            x, y = p
            return [
                (R + r*cos(x))*cos(y),
                (R + r*cos(x))*sin(y),
                r*sin(x)
            ]
        end
    end


    @with_kw struct SphereEmbedding <: AbstractEmbedding
        φ::Function = (p) -> begin
            lon, lat = p
            ls = atan(tan(lat)) 
            return [ 
                    cos(ls) * cos(lon),
                    cos(ls) * sin(lon),
                    sin(ls),
            ]
        end
    end


    @with_kw struct MobiusEmbedding <: AbstractEmbedding
        φ::Function = (p) -> begin
            x, y = p
            return [ 
                    (1 + y/2*cos(x/2))*cos(x),
                    (1 + y/2*cos(x/2))*sin(x),
                    y/2*sin(x/2)
            ]
        end
    end



    # ---------------------------------------------------------------------------- #
    #                                    VISUALS                                   #
    # ---------------------------------------------------------------------------- #

    theme = Theme(
        Axis3d = (
            backgroundcolor = :black,
            xgridcolor = :white,
            ygridcolor = :white,
        )
    )


    """
    See axis docs: https://docs.juliahub.com/MakieGallery/Ql23q/0.2.17/axis.html#Axis3D
    """

    function visualize_manifold(
        x::AbstractVector,
        y::AbstractVector,
        e::AbstractEmbedding; 
        color=nothing,
        cmap=Reverse(:bone_1),
        colorby::Symbol=:d,
        transparency::Bool=false,
    )
        X,Y,Z = e(x, y)
        
        if isnothing(color)
            if colorby == :d
                color = sqrt.(X .^ 2 .+ Y .^ 2 .+ Z .^ 2)
            elseif colorby == :z
                color = Z
            else
                error("Unrecognized colorby value $colorby")
            end
        end

        # plot
        fig = Figure(resolution=(1200, 1200), viewmode = :fitzoom)
        ax = LScene(fig[1, 1], show_axis=true, )

        pltobj = surface!(
            ax,
            X, Y, Z;
            shading=false,
            color=color,
            colormap=cmap,
            transparency=transparency,
        )
        wireframe!(ax, X, Y, Z; transparency=transparency, shading=false, color=:black, linewidth=0.5)

        # colorbar
        zoom!(ax.scene, cameracontrols(ax.scene), 1.4)
        Colorbar(fig[1, 2], pltobj, height=Relative(0.5),
            label = string(colorby), ticklabelsize = 18,
            ticklabelcolor=:white, tickcolor=:white,
            labelcolor=:white, labelsize=20,
        )
        colsize!(fig.layout, 1, Aspect(1, 0.8))
        colsize!(fig.layout, 2, Aspect(1, 0.1))


        # style
        axis = ax.scene[OldAxis] # you can change more colors here!
        axis[:ticks][:textcolor] = :grey64
        axis[:ticks][:textsize] = 12
        axis[:frame][:linecolor] = :white
        axis[:frame][:axiscolor] = :white
        axis[:frame][:linewidth] = 0.5
        axis[:names][:textcolor] = :white
        axis[:names][:textsize] = 22


        set_theme!(backgroundcolor=colorant"#23272E")
        display(fig)
    end
end


