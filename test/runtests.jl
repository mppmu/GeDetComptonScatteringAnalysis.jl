# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

using Pkg
Pkg.Registry.add(Pkg.RegistrySpec(url = "https://github.com/JuliaRegistries/General"))
Pkg.Registry.add(Pkg.RegistrySpec(url = "https://github.com/legend-exp/LegendJuliaRegistry"))
Pkg.add(url = "https://github.com/legend-exp/LegendHDF5IO.jl", rev = "main")

using Test
using GeDetComptonScatteringAnalysis

datapath = joinpath(@__DIR__, "testdata")
destdir = joinpath(@__DIR__, "results")
!isdir(destdir) && mkdir(destdir)

@testset "segBEGe" begin
    hv = 300.        # applied Voltage
    phi = 88.5       # azimuthal position of radiation source (in degrees)
    z = 62.8         # height of czt cameras (in mm)
    name = "segBEGe" # group name of detector in raw lh5 files

    @testset "only core" begin
        r = 81.8     # radial position of radiation source (in mm)
        stack_and_merge_at_z(datapath, destdir, r, phi, z, hv, name)
        resultfile = joinpath(destdir, "R_81.8mm_Z_62.8mm_Phi_88.5deg_T_71.75K_measuretime_1200sec_HV_300V-20230126T145729Z-preprocessed.lh5")
        @test isfile(resultfile)
        mtime, R, Z = get_all_z(destdir)
        rm(resultfile)
    end

    # it looks like all files exceed the 100MB/file limit, we might need to create dedicated test files for this
    #=
    @testset "core and segments" begin
        r = 79.8     # radial position of radiation source (in mm)
        stack_and_merge_at_z(datapath, destdir, r, phi, z, hv, name)
        resultfile = joinpath(destdir, "R_79.8mm_Z_62.8mm_Phi_88.5deg_T_71.75K_measuretime_600sec_HV_300V-20230118T010831Z-preprocessed.lh5")
        @test isfile(resultfile)
        mtime, R, Z = get_all_z(destdir)
        rm(resultfile)
    end
    =#
end

@testset "ICPC" begin
    # Add a test-set from ICPC here
end
