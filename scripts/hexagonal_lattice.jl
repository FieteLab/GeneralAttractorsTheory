
using Plots

"""
Playing around with hexagonal lattices
"""

w, h = 100, 100

function make_lattice(λ)::NamedTuple
    pts = (; x=[0.0,], y=[0.0,])

    p = [0, 0]
    for θ in range(0, 2π, length=7)
        θ == 0 && continue
        push!(pts.x, λ*cos(θ))
        push!(pts.y, λ*sin(θ))
    end


    x₀, y₀ = copy(pts.x), copy(pts.y)

    for Δy in (1, -1)
        append!(
            pts.x, x₀ .+ 0
        )
        append!(
            pts.y, y₀ .+ 2sin(π/3)*λ*Δy
        )
    end


    for ϕ in (0, π/3, π/3)
        Δx = cos(ϕ + π/3)

        for Δy in (1, -1)
            append!(
                pts.x, x₀ .+  Δx*λ .+ λ*(sign(Δx))
            )
            append!(
                pts.y, y₀ .+ sin(π/3)*λ*Δy
            )
        end

        
    end

    x₀, y₀ = copy(pts.x), copy(pts.y)

    append!(pts.x, x₀ .+ 4λ)
    append!(pts.y, y₀)
    append!(pts.x, x₀ .- 4λ)
    append!(pts.y, y₀)

    x₀, y₀ = copy(pts.x), copy(pts.y)
    append!(pts.x, x₀ .+ 2λ)
    append!(pts.y, y₀ .+ 6λ * sin(π/3) )
    append!(pts.x, x₀ .- 2λ)
    append!(pts.y, y₀ .- 6λ * sin(π/3))


    x₀, y₀ = copy(pts.x), copy(pts.y)
    append!(pts.x, x₀ .+ 12λ)
    append!(pts.y, y₀ )

    x₀, y₀ = copy(pts.x), copy(pts.y)
    append!(pts.x, x₀ .- 12λ)
    append!(pts.y, y₀ )

    x₀, y₀ = copy(pts.x), copy(pts.y)
    append!(pts.x, x₀)
    append!(pts.y, y₀ .+ 12λ* sin(π/3))

    x₀, y₀ = copy(pts.x), copy(pts.y)
    append!(pts.x, x₀)
    append!(pts.y, y₀ .- 12λ* sin(π/3))

    clean = zip(pts.x, pts.y) |> collect |> unique!

    pts = (;
        x = [x[1] for x in clean],
        y = [x[2] for x in clean],
    )
    return pts
end

λ = 0.2
φ1 = make_lattice(1)
φ2 = make_lattice(1 + λ*1)
φ3 = make_lattice(1 + λ*2)
φ4 = make_lattice(1 + λ*3)

scatter(φ1.x, φ1.y, ms=4, aspect_ratio=:equal, 
            size=(1000, 1000), 
            xlim=[-30, 30], ylim=[-30, 30],
        )
scatter!(φ2.x, φ2.y, ms=4)
scatter!(φ3.x, φ3.y, ms=4)
scatter!(φ4.x, φ4.y, ms=4)