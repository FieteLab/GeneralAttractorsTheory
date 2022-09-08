# GeneralAttractors

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://FedeClaudi.github.io/GeneralAttractors.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://FedeClaudi.github.io/GeneralAttractors.jl/dev/)
[![Build Status](https://github.com/FedeClaudi/GeneralAttractors.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/FedeClaudi/GeneralAttractors.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/FedeClaudi/GeneralAttractors.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/FedeClaudi/GeneralAttractors.jl)


## Folder structure

The `GeneralAttractors.jl` package code is defined in `src`.

`scripts` holds useful `.jl` files (e.g. to visualize distance metrics or pre-defined CANs).

## Usage
to define a `CAN` you need to first define the number of dimensions and the number of neurons in each dimensions. This is done with a `Tuple` $n$ of `Int` with the number f neurons in each dimensions. 

Next, you need a map `ξ` assigning to each neuron a location in the neural lattice based on it's index. 
For one-dimensionals CANs, $\xi: [1:n] \to \mathbb R$.
For two dimensional CAns: $\xi: [1:n_1]\times[1:n_2] \to \mathbb R^2$ and so on. 
The neurons positions is used to compute the pairwise distance matrix based on a metric `d`. 

Metrics are defined based on the `Distances.jl` package which includes a `PeriodicEuclidean` metric type that computes distances on a periodic domain. For the ring attractor the domain is one dimensional and with extend $2\pi$, for the torus it's $[0:2\pi] \times [0:2\pi]$ and for the cylinder it's $[0, 2\pi] \times [0, 1]$.
For the Mobius attractor a custom `MobiusEucildean` metric is defined in `src/metrics.jl`. When constructing a `CAN`, in addition to passin $n, \xi$ you also pass a metric $d$ (a subtype of `Distances.Metric`). 

Finally, the last step is to use a `Kernel` to go from distance to connectivity strength. `kernels.jl` defines a few `AbstractKernel` types, but essentially a kernel is entirely defined by a map `k: \mathbb R \to \mathbb R`, the `AbstractKernels` types are used to conveniently pass arguments to $k$ in a way that increases code performance. 


### Under the hood.
Once `CAN` is passed $n, ξ, d, K$ (in code: `Tuple, Function, Metric, AbstractKernel`) it creates a network:
- it first layous out the neurons in a lattice of shape $n$ with neurons at coordinates given by $ξ$.
- it computes the pairwise distance matrix between neurons. 
- it uses the kernel $k$ to get the connection strength.

The last two steps repeated $2d$ times with $d=length(n)$ the number of dimensions of the attractor to create $W_i$ offset connectivity matrices. To compute the offset, given tge the coordinates of all neurons in the lattice $X$ and an offset vector $Δx ∈ \mathbb R^d$, the pairwise distance matrice is computed as `metric(X + Δx, X)`.
For each matrix $\Delta x$ is choen to offset the distances by $\pm$ a "basis" of the lattice by a an offset fctor set by the `offset_strength` paramter for `CAN`. 

## Note
Currently `Kernel`s are considered to be radially symmetric, hence why they work as maps `k: \mathbb R \to \mathbb R` and why we can create the $W_i$ by offetting the pairwise distance matrix. Shuld we need asymmetric connection kernels this would have to be changed. However this solution is the most compiutationally efficient for now. 

## TODO
- [x] connect $W_i$
- [ ] simulate
- [ ] visualize activity
- [ ] PCA/isomap/TDA -> visualize and reconstruct activity manifold