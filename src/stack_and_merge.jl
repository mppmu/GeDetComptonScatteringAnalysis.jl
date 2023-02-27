# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

_vcat!(t1::Tab, t2::Tab) where {Tab <: Union{detTable, cztTable}} = begin
    t2.evt_no .+= t1.evt_no[end]
    append!(t1, t2)
end
_vcat!(t1::Tuple{detTable, cztTable}, t2::Tuple{detTable, cztTable}) = 
    (_vcat!(t1[1], t2[1]), _vcat!(t1[2], t2[2]))
_vcat!(x::Tuple{Tuple{Missing, Missing}, Int}, ::Tuple{Missing, Missing}) = x
_vcat!(::Tuple{Tuple{Missing, Missing}, Int}, y::Tuple{detTable, cztTable}) = (y, 1)
_vcat!(x::Tuple{Tuple{detTable, cztTable}, Int}, ::Tuple{Missing, Missing}) = x
_vcat!(x::Tuple{Tuple{detTable, cztTable}, Int}, y::Tuple{detTable, cztTable}) = 
    (_vcat!(x[1], y), x[2] + 1)


function read_and_merge_filtered_file(
        file::AbstractString, name::AbstractString; #hv::Number; 
        idx_c::Int = 1, #corr_daq_energy::Bool = false, rm_pileup::Bool = false
    )::Tuple{detTable, cztTable}

    det::detTable, czt1::cztTable, czt2 = read_filtered_file(file, name)
    core::detTable = det[findall(det.chid .== idx_c)]
    motor_z::QuantityMM{Float64} = getZ(file)
    czt::cztTable = merge_cameras_and_transform_coordinates(core, czt1, czt2, motor_z)
    # TODO: reintroduce rm_pileup with idx that also works for detectors with multiple channels
    # rm_pileup && (idx = intersect(findall(is_not_peak_pileup, core.samples), idx))
    det, czt
end

function stack_and_merge_at_z(
        sourcedir::AbstractString, destdir::AbstractString, 
        r::Number, phi::Number, z::Number, hv::Number, name::AbstractString; 
        idx_c::Int = 1, #corr_daq_energy::Bool = false, rm_pileup::Bool = false,
        hv_in_filename::Bool = false, n_max_files::Int = -1, verbose::Bool = true
    )::Nothing 

    files = fetch_relevant_filtered_files(sourcedir, phi, z, r, hv_in_filename, n_max_files)
    x = ((missing, missing), 0)
    for file in files
        try
            x = _vcat!(x, 
            read_and_merge_filtered_file(file, name; idx_c))
        catch e
            @warn "$file was ignored due to possible file issues"
        end
    end
    verbose && println("$(x[2]) / $(length(files)) successful")
    fileout = build_preprocessed_file_name(files, destdir, x[2])
    LHDataStore(f -> write_preprocessed_file(f, name, x[1]), fileout, "w")
    chmod(fileout, 0o774)
    nothing
end