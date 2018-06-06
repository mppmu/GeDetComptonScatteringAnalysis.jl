# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

__precompile__(true)

module GeDetComptonScatteringAnalysis

using Compat
using Compat.Markdown
using Compat: axes

using ArraysOfArrays
using DataFrames
using H3DPolaris
using IntervalSets
using LsqFit
using RecipesBase

include("multifind.jl")
include("functions.jl")
include("plot_recipes.jl")

end # module
