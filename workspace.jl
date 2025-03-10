# using Pkg
# ENV["PYTHON"] = ""  # This will use the default Python or let Julia download its own
# Pkg.build("PyCall")
# Pkg.build("PyPlot")

using Pkg
Pkg.build("CMake")
Pkg.build("MiniQhull")
Pkg.build("Ripserer")