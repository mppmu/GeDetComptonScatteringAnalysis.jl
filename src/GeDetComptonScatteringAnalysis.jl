# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

__precompile__(true)

module GeDetComptonScatteringAnalysis

using LinearAlgebra
using Statistics

using ArraysOfArrays
using ArraysOfArrays: no_consistency_checks
using DSP
using IntervalSets
using LegendHDF5IO
using LsqFit
using ProgressMeter
using RadiationDetectorDSP
using RadiationDetectorDSP: SamplesOrWaveform, RealQuantity, RDWaveform, 
    AbstractSamples, _floattype
using RadiationSpectra
using RecipesBase
using Rotations
using StaticArrays
using StatsBase
using StructArrays
using Tables
using TypedTables
using Unitful

export get_z_from_2_hit_events, get_z_from_energies, swap_CZT_hits
export polaris_dither, correct_timestamps!, find_common_events
export compton_angle, compton_E_out
export stack_and_merge_at_z, get_all_z
export get_z_and_waveforms, reconstruct_at_radius
export get_daqe
export Cone


include("types.jl")
include("functions.jl")
include("cone.jl")
include("compton_scanner_analysis.jl")
include("pileup.jl")
include("energies.jl")
include("IO.jl")
include("stack_and_merge.jl")
include("transformations.jl")
include("z_reconstruction.jl")
include("pulse_shape_analysis.jl")
# include("plot_recipes.jl")

end # module
