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

export get_z_from_2_hit_events, get_z_from_energies, swap_CZT_hits
export polaris_dither, correct_timestamps!, find_common_events
export compton_angle, compton_E_out
export stack_and_merge_at_z, get_all_z
export Cone


include("multifind.jl")
include("functions.jl")
# include("plot_recipes.jl")
include("cone.jl")
include("compton_scanner_analysis.jl")
include("pileup.jl")
include("energies.jl")
include("IO.jl")
include("stack_and_merge.jl")
include("transformations.jl")
include("new_funcs.jl")

end # module
