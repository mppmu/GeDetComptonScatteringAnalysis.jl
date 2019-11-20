# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

__precompile__(true)

module GeDetComptonScatteringAnalysis

using LinearAlgebra
using Statistics

using ArraysOfArrays
using H3DDetectorSystems
using IntervalSets
using LsqFit
using Plots # !!! Temporary
using RecipesBase
using Tables
using TypedTables
using Unitful

include("multifind.jl")
include("functions.jl")
include("plot_recipes.jl")
include("cone.jl")
include("compton_scanner_analysis.jl")

end # module
