
"""
Get the location of the peak of activity on a neuron lattice
in lattice coordinates
"""
function decode_peak_location(s::Vector, can::AbstractCAN)
    idx = argmax(s)
    return can.X[:, idx]
end


mutable struct Decoder
    x::Vector # last state in can.C.M, the cover space mfld coordinates
    n::Vector # last state in can.C.N, the neural lattice coordinates
end


"""
    (dec::Decoder)(s::Vector, can::AbstractCAN)

Decode x̂ state based on the previous state `dec.x`,
the current population activity `s` and the `can` model.

If the `can`'s `CoverSpace` has the `M` and `N` spaces
being homeomorphic then decoding is just done by looking
for the peak in neural activity `s` over the neural lattice
and get the coordinates `can.X` on the neural sheet.

If the cover space is made up of two distinct manifolds, then
the change in position over the neural lattice is used to 
"move" a stored representation of the decoded variable in 
the `M` space.

When the neural state moves around a periodic manifold `N`, 
decoded is done by looking at the pre-image of the state given the
inverse of the cover space map in `M` and finding the spots
that are closest to it. 
"""
function (dec::Decoder)(s::Vector, can::AbstractCAN)
    # get the position of activity bump in neural mfld coordinates
    peak = decode_peak_location(s, can)
    can.C.M == can.C.N && return peak

    # get Δn relative to previous bump coordinates
    Δn = peak .- dec.n

    if norm(Δn) < 5
        # for small on-mfld movement, just look at the change in coordinates
        dec.x += Δn
        dec.n = peak
    else
        # for large Δn look at the pre-image of n in the variable manifold given cover map
        candidates = can.C.ρⁱ(peak...) # shape d × N

        # get the point closest to the latest decoded location
        d = map(i -> can.C.M.metric(dec.x, candidates[:, i]'), 1:size(candidates, 2))
        selected = argmin(
            d
        )

        # @info "decoding" dec.x selected candidates[:, selected]
        # plt = scatter(eachrow(candidates)..., marker_z=d)
        # scatter!([dec.x[1]], [dec.x[2]], color=:green)
        # scatter!([dec.n[1]], [dec.n[2]], color=:blue, ms=5)
        # # scatter!([candidates[1, selected]], [candidates[2, selected]], color=:black, ms=5)

        # scatter!([peak[1]], [peak[2]], color=:red)
        # display(plt)
        
        # set it as the new "decoded" position
        dec.x = candidates[:, selected]
        
        # update stored representation of peak location on neural manifold
        dec.n = peak

        
    end
    return dec.x
end