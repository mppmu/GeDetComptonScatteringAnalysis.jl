# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

__precompile__(true)

module GeDetComptonScatteringAnalysis

using LinearAlgebra
using Statistics

using ArraysOfArrays
using IntervalSets
using LsqFit
using RecipesBase
using Tables
using TypedTables
using Unitful
using LegendHDF5IO
using Rotations
using StaticArrays
using ArraysOfArrays: no_consistency_checks


include("multifind.jl")
include("functions.jl")
# include("plot_recipes.jl")
include("cone.jl")
include("compton_scanner_analysis.jl")
include("new_funcs.jl")

end # module
