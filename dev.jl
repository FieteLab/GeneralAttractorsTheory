using GeneralAttractors

import Base.Iterators: product as ×  # cartesian product
import Distances: Metric, PeriodicEuclidean, pairwise
import StaticArrays: SVector, SA_F64, SMatrix
using Term.Progress
using Plots



using Logging
import Term.Logs: TermLogger
Logging.min_enabled_level(::TermLogger) = Logging.Debug






# ----------------------------------- ring ----------------------------------- #
@info "doing ring"


# ----------------------------------- torus ---------------------------------- #
@info "doing torus"


show_connectivity(ring_attractor,)

# TODO avoid having to bake kernel params in.


