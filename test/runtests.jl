# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

using Test
using GeDetComptonScatteringAnalysis
using GeDetComptonScatteringAnalysis: get_econv, econv
using Unitful

datapath = joinpath(@__DIR__, "testdata")
destdir = joinpath(@__DIR__, "results")
!isdir(destdir) && mkdir(destdir)
for file in readdir(destdir) rm(joinpath(destdir, file)) end

@testset "segBEGe" begin
    hv = 300.        # applied voltage (in V)
    phi = 88.5       # azimuthal position of radiation source (in degrees)
    z = 62.8         # height of czt cameras (in mm)
    name = "segBEGe" # group name of detector in raw lh5 files

    @testset "only core, czt and czt2" begin
        r = 81.8     # radial position of radiation source (in mm)
        stack_and_merge_at_z(datapath, destdir, r, phi, z, hv, name)
        resultfile = joinpath(destdir, "R_81.8mm_Z_62.8mm_Phi_88.5deg_T_71.75K_measuretime_1200sec_HV_300V-20230126T145729Z-preprocessed.lh5")
        @test isfile(resultfile)
        mtime, R, Z = get_all_z(destdir)
        @inferred reconstruct_at_radius(resultfile, hv, ew = 20u"keV")
        # get_econv is inconsistent with the values quoted in econv dictionary
        @test_broken get_econv(datapath) ≈ econv[hv]*u"keV" rtol=0.1 
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
    hv = 1200.       # applied voltage (in V)
    phi = 69.9       # azimuthal position of radiation source (in degrees)
    z = 25.0         # height of czt cameras (in mm)
    name = "ICPC"    # group name of detector in raw lh5 files

    # G.x_CZT, G.y_CZT, G.z_CZT = 1000u"μm" .* (-6.404541969299316, 70.71229063689407 + 2.5, 72.9764303564702 - 87)

    @testset "only core and czt" begin
        r = 71.0     # radial position of radiation source (in mm)
        stack_and_merge_at_z(datapath, destdir, r, phi, z, hv, name, idx_c = 16)
        resultfile = joinpath(destdir, "R_71.0mm_Z_25.0mm_Phi_69.9deg_T_95.0K_measuretime_600sec_HV_1200.0V-20220919T141651Z-preprocessed.lh5")
        @test isfile(resultfile)
        mtime, R, Z = get_all_z(destdir, name = name, center = 82.12404619872268)
        @inferred reconstruct_at_radius(resultfile, hv, name = name, idx_c = 16)
        @test get_econv(datapath, idx_c = 16, bsize = 1000, max=1_500_000, name = name) ≈ econv[hv]*u"keV" rtol=0.1
        rm(resultfile)
    end
end
