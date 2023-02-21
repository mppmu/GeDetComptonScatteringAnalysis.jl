

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


function read_and_merge_filtered_file(file, name, hv; idx_c=1, corr_daq_energy=false, rm_pileup=false)
    try
        det, czt, czt2 = read_filtered_file(file, name)
        nseg = length(unique(det.chid))
        core = det[findall(det.chid .== idx_c)]
        cf = econv[hv]
        czt = merge_data(core, czt, czt2, getZ(file))
        corr_daq_energy && correct_DAQ_energy!(core, cf)
        idx = findall(e -> 250 ≤ cf*e ≤ 440, core.DAQ_energy)
        rm_pileup && (idx = intersect(findall(is_not_peak_pileup, core.samples), idx))
        core[idx], czt[idx]
    catch  
        @warn "$file was ignored due to possible file issues"
        (missing, missing)
    end
end

function stack_and_merge_at_z(sourcedir::String, destdir::String, r, phi, 
        z, hv, name; idx_c=1, hv_in_filename=false, 
        corr_daq_energy=false, rm_pileup=false, n_max_files=-1, verbose::Bool = true)
    files = fetch_relevant_filtered_files(sourcedir, phi, z, r, hv_in_filename, 
        n_max_files)
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
end

function merge_data(::Table, czt::Table, czt2::Missing, z::Float64)
    transform_czt1_coords!(czt.hit_x.data, czt.hit_y.data, czt.hit_z.data, z)
    czt
end

function merge_data(det, czt, czt2, z)
    # TODO: maybe consider doing transformation directly in main loop here
    # instead of seperately for each camera
    transform_czt1_coords!(czt.hit_x.data, czt.hit_y.data, czt.hit_z.data, z)
    transform_czt2_coords!(czt2.hit_x.data, czt2.hit_y.data, czt2.hit_z.data, z)
    N, N1, N2 = length(det), length(czt), length(czt2)
    total_cam_hits = length(czt.hit_x.data) + length(czt2.hit_x.data)
    # simple Vectors
    evt_nhits = Vector{Int32}(undef, N)
    evt_t = Vector{eltype(czt.evt_t)}(undef, N)
    elem_ptr = Vector{Int32}(undef, N+1)
    # data for VectorOfVectors
    detno = Vector{Int32}(undef, total_cam_hits)
    edep = Vector{eltype(czt.hit_edep.data)}(undef, total_cam_hits)
    t = Vector{eltype(czt.hit_t.data)}(undef, total_cam_hits)
    x = Vector{eltype(czt.hit_x.data)}(undef, total_cam_hits)
    y = Vector{eltype(czt.hit_y.data)}(undef, total_cam_hits)
    z = Vector{eltype(czt.hit_z.data)}(undef, total_cam_hits)

    i1, i2, ep = 1, 1, 0
    for i=eachindex(det)
        elem_ptr[i] = ep + 1
        evt_no_cam1, evt_no_cam2 = czt.evt_no[i1], czt2.evt_no[i2]
        # events with hits in both cameras
        if evt_no_cam1 == evt_no_cam2 == det.evt_no[i]
            evt_nhits[i] = czt.evt_nhits[i1] + czt2.evt_nhits[i2]
            evt_t[i] = min(czt.evt_t[i1], czt2.evt_t[i2])
            for n=Base.OneTo(czt.evt_nhits[i1])
                detno[ep + n] = czt.hit_detno[i1][n]
                edep[ep + n] = czt.hit_edep[i1][n]
                t[ep + n] = czt.hit_t[i1][n]
                x[ep + n] = czt.hit_x[i1][n]
                y[ep + n] = czt.hit_y[i1][n]
                z[ep + n] = czt.hit_z[i1][n]
            end
            ep += czt.evt_nhits[i1]
            for n=Base.OneTo(czt2.evt_nhits[i2])
                detno[ep + n] = 4 + czt2.hit_detno[i2][n]
                edep[ep + n] = czt2.hit_edep[i2][n]
                t[ep + n] = czt2.hit_t[i2][n]
                x[ep + n] = czt2.hit_x[i2][n]
                y[ep + n] = czt2.hit_y[i2][n]
                z[ep + n] = czt2.hit_z[i2][n]
            end
            ep += czt2.evt_nhits[i2]
            i1 += 1
            i2 += 1
        # events with hits in only camera 1
        elseif evt_no_cam1 == det.evt_no[i]
            evt_nhits[i] = czt.evt_nhits[i1]
            evt_t[i] = czt.evt_t[i1]
            for n=Base.OneTo(czt.evt_nhits[i1])
                detno[ep + n] = czt.hit_detno[i1][n]
                edep[ep + n] = czt.hit_edep[i1][n]
                t[ep + n] = czt.hit_t[i1][n]
                x[ep + n] = czt.hit_x[i1][n]
                y[ep + n] = czt.hit_y[i1][n]
                z[ep + n] = czt.hit_z[i1][n]
            end
            ep += czt.evt_nhits[i1]
            i1 += 1
        # events with hits in only camera 2
        elseif evt_no_cam2 == det.evt_no[i]
            evt_nhits[i] = czt2.evt_nhits[i2]
            evt_t[i] = czt2.evt_t[i2]
            for n=Base.OneTo(czt2.evt_nhits[i2])
                detno[ep + n] = 4 + czt2.hit_detno[i2][n]
                edep[ep + n] = czt2.hit_edep[i2][n]
                t[ep + n] = czt2.hit_t[i2][n]
                x[ep + n] = czt2.hit_x[i2][n]
                y[ep + n] = czt2.hit_y[i2][n]
                z[ep + n] = czt2.hit_z[i2][n]
            end
            ep += czt2.evt_nhits[i2]
            i2 += 1
        # events with no hits in the cameras
        # this should not happen, as the filtered files should
        # only contain events with hits in the detector AND the camera
        else 
            @warn "Event with event number $(det.evt_no[i]) in the detector has no counterpart in the cameras."
        end
        i1 > N1 && (i1 = N1)            # TODO: maybe get rid of these?
        i2 > N2 && (i2 = N2)
    end
    elem_ptr[end] = total_cam_hits + 1
    Table(
        evt_no = det.evt_no,
        evt_nhits = evt_nhits,
        evt_t = evt_t,
        hit_detno = VectorOfVectors(detno, elem_ptr, no_consistency_checks),
        hit_edep = VectorOfVectors(edep, elem_ptr, no_consistency_checks),
        hit_t = VectorOfVectors(t, elem_ptr, no_consistency_checks),
        hit_x = VectorOfVectors(x, elem_ptr, no_consistency_checks),
        hit_y = VectorOfVectors(y, elem_ptr, no_consistency_checks),
        hit_z = VectorOfVectors(z, elem_ptr, no_consistency_checks),
        )
end