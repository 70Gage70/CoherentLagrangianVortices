import Pkg
Pkg.activate(@__DIR__)

Pkg.add(url = "https://github.com/CoherentStructures/CoherentStructures.jl")
Pkg.add(["ThreadsX", "MAT", "StaticArrays", "AxisArrays", 
    "Dates", "Plots", "NetCDF", "NCDatasets", "Interpolations", 
    "ScatteredInterpolation", "OrdinaryDiffEq", "Tensors", "MKL_jll"])

Pkg.update()