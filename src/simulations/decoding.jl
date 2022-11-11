
"""
Get the location of the peak of activity on a neuron lattice
in lattice coordinates
"""
function decode_peak_location(s::Vector, can)
    idx = argmax(s)
    return can.X[:, idx]
end