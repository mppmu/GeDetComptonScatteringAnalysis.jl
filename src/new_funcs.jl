export stack_and_merge_at_z, get_all_z

regR = "R_(?<R>\\d+?\\.\\d+)mm"
regPhi = "_Phi_(?<Phi>\\d+?\\.\\d+)deg"
regT = "_T_(?<T>\\d+?\\.\\d+)K"
regZ = "_Z_(?<Z>\\d+?\\.\\d+)mm"
regM = "_measuretime_(?<time>\\d+?)sec"
regV = "HV_(?<voltage>\\d+?)V"
getR(s::String) = parse(Float64, match(Regex(regR), s).captures[1])
getPhi(s::String) = parse(Float64, match(Regex(regPhi), s).captures[1])
getT(s::String) = parse(Float64, match(Regex(regT), s).captures[1])
getZ(s::String) = parse(Float64, match(Regex(regZ), s).captures[1])
getM(s::String) = parse(Float64, match(Regex(regM), s).captures[1])
getV(s::String) = parse(Float64, match(Regex(regV), s).captures[1])

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

const econv = Dict( 
                    300.0 => 0.000746364303766692,
                    # 600.0 => 0.0007730410437211994,
                    600.0 => 0.0006952481325306902,
                    # 900.0 => 0.0006666592487136949,
                    900.0 => 0.0006781005477809491,
                    1200.0 => 0.0006189596091123743, 
                    3500.0 => 0.0005490928284752016,
                    3000.0 => 0.000653
                )
const x_CZT, y_CZT, z_CZT = 1000u"μm" .* (-4.5, 81.764-1.55, -62.8 + 19.75)
const δ = [x_CZT, y_CZT, z_CZT]
const α2 = 45.3695 * π/180
const δ2 = 1000u"μm" .* [64.1275058422767, -31.859041946909453, 0.4]

write_all(f::LHDataStore, name::String , data::Tuple{Table, Table}) = begin
    LegendHDF5IO.writedata(f.data_store, "$name", data[1])
    LegendHDF5IO.writedata(f.data_store, "czt", data[2])
end

function stack_and_merge_at_z(sourcedir::String, destdir::String, r, phi, 
        z, hv, name; idx_c=1, hv_in_filename=false, 
        corr_daq_energy=false, rm_pileup=false, n_max_files=-1)
    files = fetch_relevant_files(sourcedir, phi, z, r, hv_in_filename, 
        n_max_files)
    _read_file(file) = read_file(
        file, name, hv; idx_c, corr_daq_energy, rm_pileup)
    x = ((missing, missing), 0)
    for i=eachindex(files)
        println("at $i")
        x = _vcat!(x, _read_file(files[i]))
    end
    println("$(x[2]) / $(length(files)) successfull")
    fileout = build_output_fname(files, destdir, x[2])
    LHDataStore(f -> write_all(f, name, x[1]), fileout, "w")
    chmod(fileout, 0o774)
end

function fetch_relevant_files(sourcedir, phi, z, r, hv_in_fname, n_max_files)
    sourcedir = joinpath(sourcedir, "")
    allfiles = filter(file -> endswith(file, "filtered.h5")
                            && occursin("Phi_$(phi)", file)
                            && occursin("Z_$(z)", file)
                            && occursin("R_$(r)", file), readdir(sourcedir))
    if hv_in_fname
        allfiles = filter(file -> contains(file, "_HV_$hv"), allfiles)
    end
    if n_max_files > 0
        allfiles = allfiles[1:min(length(allfiles), n_max_files)]
    end
    joinpath.(sourcedir, allfiles)
end

function build_output_fname(files, destdir, n)
    destdir = joinpath(destdir, "")
    fpos, ftime = split(basename(files[end]), "measuretime_")
    fmtime, fdate = split(ftime, "sec")
    fmtime = n*parse(Int32, fmtime)
    joinpath(destdir, fpos * "measuretime_$(fmtime)sec" * fdate)
end

function get_nseg(f::LHDataStore, name::String)
    chid = f[name*"/chid"][1:10]
    length(unique(chid))
end

function read_file(file, name, hv; idx_c=1, corr_daq_energy=false, 
        rm_pileup=false)
    try
        _get_nseg(lhd) = get_nseg(lhd, name)    # TODO: should we stay flexible?
        _read_all(lhd) = read_all(lhd, name)
        nseg = LHDataStore(_get_nseg, file)
        detector, czt, czt2 = LHDataStore(_read_all, file)
        core = detector[idx_c:nseg:length(detector)]
        cf = econv[hv]
        czt = merge_data(core, czt, czt2, getZ(file))
        corr_daq_energy && correct_DAQ_energy!(core, cf)
        idx = findall(e -> 250 ≤ cf*e ≤ 440, core.DAQ_energy)
        rm_pileup && 
            (idx = intersect(findall(is_not_peak_pileup, core.samples), idx))
        core[idx], czt[idx]
    catch  
        @warn "$file was ignored due to possible file issues"
        (missing, missing)
    end
end

read_all(f::LHDataStore, name::String) = 
    (f[name][:], f["czt"][:], "czt2" in keys(f) ? f["czt2"][:] : missing)

function merge_data(::Table, czt::Table, czt2::Missing, z::Float64)
    transform_czt1_coords!(czt.hit_x, czt.hit_y, czt.hit_z, z)
end

function merge_data(ge, czt, czt2, z)
    # TODO: maybe consider doing transformation directly in main loop here
    # instead of seperately for each camera
    transform_czt1_coords!(czt.hit_x.data, czt.hit_y.data, czt.hit_z.data, z)
    transform_czt2_coords!(czt2.hit_x.data, czt2.hit_y.data, czt2.hit_z.data, z)
    N, N1, N2 = length(ge), length(czt), length(czt2)
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
    for i=eachindex(ge)
        elem_ptr[i] = ep + 1
        evt_no_cam1, evt_no_cam2 = czt.evt_no[i1], czt2.evt_no[i2]
        if evt_no_cam1 == evt_no_cam2 == ge.evt_no[i]
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
        elseif  evt_no_cam1 == ge.evt_no[i]
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
        else
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
        end
        i1 > N1 && (i1 = N1)            # TODO: maybe get rid of these?
        i2 > N2 && (i2 = N2)
    end
    elem_ptr[end] = total_cam_hits + 1
    Table(
        evt_no = ge.evt_no,
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

# lowest allocation method
function transform_czt2_coords!(X::Vector{T}, Y::Vector{T}, Z::Vector{T}, 
motor_z::U) where {T, U}
    t = SVector{3}(1000u"μm" .* (U(0), U(0), motor_z) .+ δ .+ δ2)
    R = RotZ(-α2)*RotY(π)
    @inbounds for i=eachindex(X)
        v = SVector{3, eltype(t)}(X[i], -Z[i], -Y[i])
        X[i], Y[i], Z[i] = round.(T, R * v + t)
    end
end

# lowest allocation method 
function transform_czt1_coords!(X::Vector{T}, Y::Vector{T}, Z::Vector{T}, 
motor_z::U) where {T, U}
    t = round.(T, (1000u"μm" .* (U(0), U(0), motor_z) .+ δ))
    @inbounds for i=eachindex(X)
        Y_i = Y[i]
        X[i] = t[1] + X[i]
        Y[i] = t[2] - Z[i]
        Z[i] = t[3] - Y_i
    end
end

# ----- INFO ------
# Energies should be determined from pulses...
function correct_DAQ_energy!(det, cf)
    for i=eachindex(det)
        new_e = DAQ_energy_corr(det.samples[i], cf*det.DAQ_energy[i])
        det.DAQ_energy[i] = new_e / cf
    end
end

function DAQ_energy_corr(wf, daqe)
    fit_par = [0.3762989397047389, 0.7729112545305099, 0.0007296045244350602]
    bs = baseline_slope(wf)
    if bs < -1.5
        return daqe - par(bs, fit_par)
    else
        return daqe
    end
end   

@. par(x, p) = p[1] + p[2]x + p[3]x^2

baseline_slope(wf) = mean(view(wf, 1400:1500)) - mean(view(wf, 1:100))

function is_not_pileup(wf, cut)
    try
        wf = rm_baseline(wf)
        -cut < mean(view(wf, 1800:2000)) < cut && count_peaks(wf) ≤ 1 
    catch
        false
    end 
end

function is_not_peak_pileup(wf)
    wf = rm_baseline(wf)
    try
        count_peaks(wf) ≤ 1 
    catch
        false
    end 
end

function count_peaks(wf)
    xings_i = []
    xings = 0
    twf = maw(trap_filter(wf,18,250,55.0),100)
    thold = mean(view(twf,1:200)) + 3e6
    for i in 1:length(twf)-1
        if twf[i] ≤ thold && twf[i+1] ≥ thold
            push!(xings_i, i)
        end
    end
    for i in 1:length(xings_i)-1
        if xings_i[i+1] - xings_i[i] > 100
            xings += 1
        end
    end
    xings +=1 #final xing
end

function rm_baseline(wf)
    av = sum(view(wf,1:200))/200
    return wf .- av
end

function getz(file; name="segBEGe", center=81.76361317572471, ew=8)
    icpc, czt = LHDataStore(file) do lhd
        lhd[name][:], lhd["czt"][:]
    end
    hv = getV(file)
    ec = econv[hv]
    icpc_e = ec*icpc.DAQ_energy
    czt_e = ustrip(sum.(czt.hit_edep)) / 1000
    idx = intersect(findall(x -> abs(x - 662) ≤ ew, icpc_e+czt_e), 
                    findall(x -> 250 ≤ x ≤ 440, icpc_e))
    icpc_hits = view(icpc, idx)
    czt_hits = view(czt, idx)
    R = center - getR(file)

    @info "Reconstructing Z from two hit events at R = $R"
    idx_2h = findall(is_valid_2hit, czt_hits);
    czt_2hit = view(czt_hits, idx_2h);
    icpc_2hit = view(icpc_hits, idx_2h);
    # TODO: resolve allocation issues by passing fixed empty array
    zrec2hit = get_z_from_2_hit_events.(icpc_2hit, czt_2hit, R; Δz = 2, hv);
    idx_val_1 = findall(x -> x[1] == 1, zrec2hit);
    idx_val_2 = findall(x -> x[1] == 2, zrec2hit);
    idx_val = vcat(idx_val_1, idx_val_2);

    # TODO: decide on what reconstructed z we really want 
    # from core -> 2
    # from czt  -> 3
    vcat([[x[2] x[3]] for x in view(zrec2hit, idx_val)]...)
end

@inline function is_valid_2hit(evt)
    length(evt.hit_x) == 2 && hypot(evt.hit_x[2] - evt.hit_x[1],
        evt.hit_y[2] - evt.hit_y[1],
        evt.hit_z[2] - evt.hit_z[1]) > 3.0u"mm"
end

"""
    get_all_z(sourcedir; center=CNTR_VALUE)

loop through all files in `sourcedir` and return for each file the radius, 
measuretime and reconstructed z's
"""
function get_all_z(sourcedir; name="segBEGe", center=81.76361317572471)
    ffiles = readdir(sourcedir)
    R = zeros(length(ffiles))
    mtime = zeros(length(ffiles))
    z = Vector{Float64}[]
    for i=eachindex(ffiles)
        R[i] = getR(ffiles[i])
        mtime[i] = getM(ffiles[i])
        rec_zs = getz(joinpath(sourcedir, ffiles[i]); name, center, ew=20)
        push!(z, rec_zs[:, 1])
    end
    mtime, R, z
end

