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
        r::Number, phi::Number, z::Number, hv::Number, 
        det_name::AbstractString; idx_c::Int = 1, 
        #corr_daq_energy::Bool = false, rm_pileup::Bool = false,
        bsize::Int = 1000, max::Int = 1_500_000, hv_in_filename::Bool = false, 
        n_max_files::Int = -1, verbose::Bool = true)::Nothing 
    files = fetch_relevant_filtered_files(
        sourcedir, phi, z, r, hv_in_filename, n_max_files)
    det::detTable, czt::cztTable = _detTable(), _cztTable()
    successful_files = String[]
    for i=eachindex(files)
        try
            _det::detTable, _czt::cztTable = 
                read_and_merge_filtered_file(files[i], det_name; idx_c)
            if i > 1
                last_evt_no = det.evt_no[end]
                _det.evt_no .+= last_evt_no
                _czt.evt_no .+= last_evt_no
            end
            append!(det, _det)
            append!(czt, _czt)
            append!(successful_files, [files[i]])
        catch e
            nothing
        end
    end
    econv::typeof(Cs_energy) = get_econv(det; idx_c, bsize, max, verbose)
    verbose && println(
        "$(length(successful_files)) / $(length(files)) successful")
    fileout = build_preprocessed_file_name(successful_files, destdir)
    write_preprocessed_file(fileout, det_name, (det, czt, idx_c, econv))
    chmod(fileout, 0o754)
    nothing
end