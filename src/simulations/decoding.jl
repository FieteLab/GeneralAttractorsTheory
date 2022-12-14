
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
    x̂::Vector # last state in a space homeomorphic to can.C.M but with no scaling compared to can.C.N
    n::Vector # last state in can.C.N, the neural lattice coordinates
    Δ::Vector # shift between x,n at time t=0 (they will be offset when decoding start)
    α::Float64  # scaling correction factor
end

Decoder(x::Vector, n::Vector, α) = Decoder(x, x, n, x - n, α)
Decoder(x::Vector, n::Vector) = Decoder(x, x, n, x - n, 1)



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
    n̂ = decode_peak_location(s, can)
    can.C.M == can.C.N && return n̂

    # get Δn relative to previous bump coordinates
    Δn = n̂ .- dec.n

    if norm(Δn) < 1
        # for small on-mfld movement, just look at the change in coordinates
        dec.x += Δn * dec.α
        dec.x̂ += Δn
        dec.n = n̂
    else
        # for large Δn look at the pre-image of n in the variable manifold given cover map
        candidates = can.C.ρⁱ(n̂...) # shape d × N

        # get the point closest to the latest decoded location (minus manifold shift)
        d = map(
            i -> can.C.M.metric(dec.x̂ - dec.Δ, candidates[:, i]),
            1:size(candidates, 2),
        )
        selected = argmin(d)
        @debug "decoding recalibration triggered" d[selected]

        # find potential errors in decoding
        if d[selected] > 2.5
            r(x) = round(x; digits = 2)
            @warn "Decoding problems" r.(dec.x) r.(dec.n) r.(n̂) dec.Δ
        end

        # get change in latent representation with no scaling
        x̂ = candidates[:, selected] + dec.Δ
        Δx̂ = x̂ .- dec.x̂

        # set it as the new "decoded" position using scaling
        dec.x += Δx̂ * dec.α
        dec.x̂ = x̂

        # update stored representation of n̂ location on neural manifold
        dec.n = n̂
    end
    return dec.x, n̂
end
