# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).


function read_and_merge_filtered_file(file::AbstractString, 
        name::AbstractString; idx_c::Int = 1)::Tuple{detTable, cztTable}
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
        bsize::Int = 1000, max::Int = 1_500_000, hv_in_filename::Bool = false, 
        n_max_files::Int = -1, verbose::Bool = true)::Nothing 
    files = fetch_relevant_filtered_files(sourcedir, phi, z, r, hv_in_filename, n_max_files)
    det::detTable, czt::cztTable = 
        read_and_merge_filtered_file(files[1], name; idx_c)
    successful = 1
    for file in files[2:end]
        _det::detTable, _czt::cztTable = 
            read_and_merge_filtered_file(file, name; idx_c)
        append!(det, _det)
        append!(czt, _czt)
        successful += 1
    end
    econv::typeof(Cs_energy) = get_econv(det; idx_c, bsize, max, verbose)
    verbose && println("$successful / $(length(files)) successful")
    fileout = build_preprocessed_file_name(files, destdir, successful)
    LHDataStore(
        f -> write_preprocessed_file(f, name, (det, czt, idx_c, econv)), 
        fileout, "w")
    chmod(fileout, 0o754)
    nothing
end