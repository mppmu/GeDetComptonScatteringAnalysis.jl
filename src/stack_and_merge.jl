# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).


function read_and_merge_filtered_file(file::AbstractString, 
        name::AbstractString; idx_c::Int = 1)::Tuple{detTable, cztTable}
    det, czt1, czt2 = try
        read_filtered_file(file, name)
    catch e
        @warn "Error while reading the file $file"
        throw(e)
    end
    core = @view det[findall(det.chid .== idx_c)]
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
sourcedir::AbstractString, destdir::AbstractString, r::Number, phi::Number, 
z::Number, hv::Number, det_name::AbstractString; idx_c::Int = 1, 
bsize::Int = 1000, max::Int = 1_500_000, hv_in_filename::Bool = false, 
n_max_files::Int = -1, verbose::Bool = true, overwrite::Bool = false
#corr_daq_energy::Bool = false, rm_pileup::Bool = false,
)::Nothing 
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
    core = view(det, findall(det.chid .== idx_c))
    econv::typeof(Cs_energy) = get_econv(core; bsize, max, verbose)
    verbose && println(
        "$(length(successful_files)) / $(length(files)) successful")
    fileout = build_preprocessed_file_name(successful_files, destdir)
    files = readdir(destdir)
    idx = findfirst(isequal(fileout), files)

    # check if filename already exists
    if isnothing(idx)   # file does not exist
        part1, part2 = split(fileout, "measuretime_")
        _, part2 = split(part2, "sec")
        reg = Regex("(?<=$(part1*"measuretime_"))\
            \\d+(?:\\.\\d+){0,1}s(?=$("ec"*part2))")
        idx = findfirst(f -> !isnothing(match(reg, f)), files)

        # check if only measuretime changed
        if isnothing(idx)
            # no similar file exists, so just write out the file
            write_preprocessed_file(
                fileout, det_name, (det, czt, idx_c, econv))
            chmod(fileout, 0o754)
        else
            # there exists a similar file with different measuretime
            @warn "found file with similar parameters but different \
                measuretime. Old data in $(files[idx]) will be replaced"
            write_preprocessed_file(
                fileout, det_name, (det, czt, idx_c, econv), "w")
        end
    else    # file does exist
        if overwrite
            write_preprocessed_file(
                fileout, det_name, (det, czt, idx_c, econv), "w")
        else
            # append the data?
            # what about idx_c and econv
            write_preprocessed_file(
                fileout, det_name, (det, czt, idx_c, econv), "w")
        end
    end
    nothing
end