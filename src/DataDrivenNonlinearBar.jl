module DataDrivenNonlinearBar

using LinearAlgebra, SparseArrays, StaticArrays, GaussQuadrature

include("base.jl")
include("dataset.jl")
include("assembly.jl")
include("solver.jl")
# include("newtonraphsonstep.jl")


end # module DataDrivenNonlinearBar
