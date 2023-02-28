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
    det, czt1, czt2 = try
        read_filtered_file(file, name)
    catch e
        @warn "Error while reading the file $file"
        throw(e)
    end
    core::detTable = det[findall(det.chid .== idx_c)]
    motor_z::QuantityMM{Float64} = getZ(file)
    czt::cztTable = try
        merge_cameras_and_transform_coordinates(core, czt1, czt2, motor_z)
    catch e
        @warn "Error while merging czt Tables"
        throw(e)
    end
    # TODO: reintroduce rm_pileup with idx that also works for detectors with multiple channels
    # rm_pileup && (idx = intersect(findall(is_not_peak_pileup, core.samples), idx))
    det, czt
end

function stack_and_merge_at_z(
        sourcedir::AbstractString, destdir::AbstractString, 
        r::Number, phi::Number, z::Number, hv::Number, name::AbstractString; 
        idx_c::Int = 1, #corr_daq_energy::Bool = false, rm_pileup::Bool = false,
        bsize::Int = 1000, max::Int = 1_500_000,
        hv_in_filename::Bool = false, n_max_files::Int = -1, verbose::Bool = true
    )::Nothing 

    files = fetch_relevant_filtered_files(sourcedir, phi, z, r, hv_in_filename, n_max_files)
    x = ((missing, missing), 0)
    for file in files
        try
            x = _vcat!(x, 
            read_and_merge_filtered_file(file, name; idx_c))
        catch e
            @warn "$file was ignored due to error while reading the file
                or merging the czt cameras"
        end
    end
    econv::typeof(Cs_energy) = get_econv(x[1][1]; idx_c, bsize, max, verbose)
    verbose && println("$(x[2]) / $(length(files)) successful")
    fileout = build_preprocessed_file_name(files, destdir, x[2])
    LHDataStore(f -> write_preprocessed_file(f, name, (x[1]..., idx_c, econv)), fileout, "w")
    chmod(fileout, 0o754)
    nothing
end