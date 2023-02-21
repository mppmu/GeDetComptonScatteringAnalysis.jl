

_vcat!(t1::Table, t2::Table) = begin
    t2.evt_no .+= t1.evt_no[end]
    append!(t1, t2)
end
_vcat!(t1::Tuple{Table, Table}, t2::Tuple{Table, Table}) = 
    (_vcat!(t1[1], t2[1]), _vcat!(t1[2], t2[2]))
_vcat!(x::Tuple{Tuple{Missing, Missing}, Int}, ::Tuple{Missing, Missing}) = x
_vcat!(::Tuple{Tuple{Missing, Missing}, Int}, y::Tuple{Table, Table}) = (y, 1)
_vcat!(x::Tuple{Tuple{Table, Table}, Int}, ::Tuple{Missing, Missing}) = x
_vcat!(x::Tuple{Tuple{Table, Table}, Int}, y::Tuple{Table, Table}) = 
    (_vcat!(x[1], y), x[2] + 1)


function read_and_merge_filtered_file(
        file::AbstractString, name::AbstractString, hv::Number; 
        idx_c::Int = 1, corr_daq_energy::Bool = false, rm_pileup::Bool = false
    )
    try
        det, czt, czt2 = read_filtered_file(file, name)
        nseg = length(unique(det.chid))
        core = det[findall(det.chid .== idx_c)]
        cf = econv[hv]
        czt = merge_cameras_and_transform_coordinates(core, czt, czt2, getZ(file))
        corr_daq_energy && correct_DAQ_energy!(core, cf)
        idx = findall(e -> 250 ≤ cf*e ≤ 440, core.DAQ_energy)
        rm_pileup && (idx = intersect(findall(is_not_peak_pileup, core.samples), idx))
        core[idx], czt[idx]
    catch  
        @warn "$file was ignored due to possible file issues"
        (missing, missing)
    end
end

function stack_and_merge_at_z(
        sourcedir::AbstractString, destdir::AbstractString, 
        r::Number, phi::Number, z::Number, hv::Number, name::AbstractString; 
        idx_c::Int = 1, corr_daq_energy::Bool = false, rm_pileup::Bool = false,
        hv_in_filename::Bool = false, n_max_files::Int = -1, verbose::Bool = true
    )::Nothing 

    files = fetch_relevant_filtered_files(sourcedir, phi, z, r, hv_in_filename, n_max_files)
    x = ((missing, missing), 0)
    for file in files
        x = _vcat!(x, 
        read_and_merge_filtered_file(file, 
            name, hv; idx_c, corr_daq_energy, rm_pileup))
    end
    verbose && println("$(x[2]) / $(length(files)) successful")
    fileout = build_output_prepocessed_file_name(files, destdir, x[2])
    LHDataStore(f -> write_preprocessed_file(f, name, x[1]), fileout, "w")
    chmod(fileout, 0o774)
    nothing
end