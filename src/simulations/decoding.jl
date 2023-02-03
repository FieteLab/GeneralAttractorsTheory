
"""
Get the location of the peak of activity on a neuron lattice
in lattice coordinates
"""
function decode_peak_location(s, can::AbstractCAN)
    idx = argmax(s)
    return can.X[:, idx]
end


mutable struct Decoder
    x::Vector # last state in can.C.M, the cover space mfld coordinates
    n::Vector # last state in can.C.N, the neural lattice coordinates
    decoding_offset
end

Decoder(x::Vector, n::Vector; decoding_offset=zeros(length(X))) = Decoder(x, n, decoding_offset)




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
function (dec::Decoder)(s::Vector, can::AbstractCAN)::Tuple{Vector, Vector}
    # get the position of activity bump in neural mfld coordinates
    n̂ = decode_peak_location(s, can) 
    can.C.M == can.C.N && return (n̂ .- dec.decoding_offset, n̂)

    # get Δn relative to previous bump coordinates
    # Δn = can.C.ρⁱ(n̂ .- dec.n)
    Δn = n̂ .- dec.n

    # scale Δn by the cover space scaling functions
    for (i, λ) in enumerate(can.C.λs)
        Δn[i] = λ(Δn[i])
    end

    # correct for large shifts due to periodicity in neural manifold
    if norm(Δn .* can.C.N.periodic_dimensions) > 1
        needs_correction = abs.(Δn) .> 1
        Δn = Δn - sign.(Δn)  .* 
                    can.C.N.xmax .* 
                    can.C.N.periodic_dimensions .*
                    needs_correction
    end

    dec.x = apply_boundary_conditions!(dec.x + Δn, can.C.M)[1]
    dec.n = n̂
    return dec.x, n̂
end
