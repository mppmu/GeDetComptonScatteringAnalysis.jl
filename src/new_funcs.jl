export stack_and_merge_at_z, get_all_z



const x_CZT, y_CZT, z_CZT = 1000u"μm" .* (-4.5, 81.764-1.55, -62.8 + 19.75)
const δ = [x_CZT, y_CZT, z_CZT]
const α2 = 45.3695 * π/180
const δ2 = 1000u"μm" .* [64.1275058422767, -31.859041946909453, 0.4]


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

