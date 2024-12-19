import Pkg
Pkg.activate(@__DIR__)

Pkg.add(url = "https://github.com/CoherentStructures/CoherentStructures.jl")
Pkg.add(["Distributed", "MAT", "StaticArrays", "AxisArrays", 
    "Dates", "Plots", "NetCDF", "NCDatasets", "Interpolations", 
    "ScatteredInterpolation", "OrdinaryDiffEq", "Tensors"])

Pkg.update()