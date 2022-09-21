"""
Collection of code useful to visualize manifolds
"""
module ManifoldUtils
    using GLMakie
    using Parameters

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
            ls = atan(tan(lat))    # lambda
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

    function visualize_manifold(
        x::AbstractVector,
        y::AbstractVector,
        e::AbstractEmbedding; color=nothing,
        cmap=Reverse(:bone_1),
        colorby::Symbol=:d
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

        fig = Figure(resolution=(1200, 1200), fontsize=22)
        ax = LScene(fig[1, 1], show_axis=true)
        pltobj = surface!(
            ax,
            X, Y, Z;
            shading=true,
            ambient=Vec3f(0.65, 0.65, 0.65),
            backlight=1.0f0,
            color=color,
            colormap=cmap,
            transparency=false,
        )
        wireframe!(ax, X, Y, Z; transparency=false, color=:gray, linewidth=0.5)
        # cam = cameracontrols(ax.scene)
        # cam.lookat[] = [0, 0, 200] ./ 1000
        # cam.eyeposition[] = [5000, 2000, 2000] ./ 1000
        # cam.upvector[] = [0, 1, 0]
        # update_cam!(ax.scene, cam)
        # zoom!(ax.scene, cameracontrols(ax.scene), 1.1)
        display(fig)
    end
end


