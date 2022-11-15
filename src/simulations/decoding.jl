
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
    peak = decode_peak_location(s, can)
    can.C.M == can.C.N && return peak

    shift = peak .- dec.n

    if norm(shift) < 5
        dec.x += shift
        dec.n = peak
    else
        # ? triggered by a large shift -> wrapping around a periodic manifold
        candidates = can.C.ρⁱ(dec.x...)

        @info "decoding got candidates" dec.x peak shift candidates
    end
    return dec.x
end