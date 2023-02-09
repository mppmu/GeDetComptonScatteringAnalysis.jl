in_mm(x::Number) = ustrip(uconvert(u"mm", x))

function get_z_from_2_hit_events(s::NamedTuple, c::NamedTuple, R::Number; Δz = 2, hv)

    #if econv*s.DAQ_energy < 250
        #return (false, -1)
    #end
    #get z coordinates for scatter point in segBEGe
    zθ = get_z_from_energies(s, c, R, hv)
    zα = get_z_from_camera(c, R)

    #compare them: if they agree within Δz, keep them
    for z in zα 
        if abs(z - zθ) < Δz 
            return (1, zθ, z)
        end
    end
    
    #swaphits and try again
    cs = swap_CZT_hits(c)
    zθ = get_z_from_energies(s, cs, R, hv)
    zα = get_z_from_camera(cs, R)
    for z in zα 
        if abs(z - zθ) < Δz
            return (2, zθ, z)
        end
    end

    #else discard them
    return (false, zθ)
end


function get_z_from_energies(s::NamedTuple, c::NamedTuple, R,hv)
    cf = econv[hv]
    T = typeof(1.0*unit(eltype(c.hit_x)))
    # TODO: add flag if we want to rely more on camera or DAQ_energies?
    # if ge not depleted -> worse energy resolution -> rely more on camera 
    # energy resolution
    # DAQ_energies
    # θ = compton_angle(cf*s.DAQ_energy*u"keV"+sum(c.hit_edep), sum(c.hit_edep))
    # CZT energies
    θ = compton_angle(Cs_energy, sum(c.hit_edep)) 
    in_mm(T(c.hit_z[1]) + hypot(T(c.hit_x[1]), T(c.hit_y[1]) - R*u"mm") * cot(θ))
end


function get_z_from_camera(c::NamedTuple, R)
    x_global, y_global, z_global = get_global_cam_positions(c)
    α = compton_angle(sum(c.hit_edep), c.hit_edep[2])
    zα = []

    if !(isnan(α))
        try
            #define the cone and get intersections z1,z2 with the beam axis
            camhit1 = [x_global.cam[1], y_global.cam[1], z_global.cam[1]]
            camhit2 = [x_global.cam[2], y_global.cam[2], z_global.cam[2]]
            cone = Cone(in_mm.(camhit1), in_mm.(camhit1-camhit2), α)
            z1, z2 = get_possible_z_from_camera(cone, R)

            #information about sign of cos(α) is lost when calculating zα
            #so check whether the actual vectors return α (keep) or π-α (discard)
            if validate_z(z1,cone,R)
                push!(zα,z1)
            elseif validate_z(z2,cone,R)
                push!(zα,z2)
            end
        catch e
            #this catches all imaginary solutions in get_possible_z_from_camera
            if !(e isa DomainError) error(e) end
        end
    end
    return zα
end


function get_possible_z_from_camera(cone::Cone, R::Number)

    #get relevant cone parameters
    H = cone.axis
    c0 = cone.origin - Vector([0, R, 0])
    cosα = cos(cone.α)

    #solve quadartic equation for z
    A = H[3]^2 - cosα^2
    B = -2*H[3]*dot(H,c0) + 2*c0[3]*cosα^2
    C = dot(H,c0)^2 - norm(c0)^2 * cosα^2
    z1 = (-B - sqrt(B^2 - 4*A*C)) / (2*A)
    z2 = (-B + sqrt(B^2 - 4*A*C)) / (2*A)

    return z1, z2

end


function validate_z(z::AbstractFloat, cone::Cone, R::AbstractFloat; Δα::Number = 1e-6)

    #calculate the angle from the vectors
    tmp = [0,R,z] - cone.origin
    αnew = acos(dot(cone.axis, tmp/norm(tmp)))

    #keep only true angles
    if abs(cone.α - αnew) < Δα
        return true
    else
        return false
    end

end

function swap_CZT_hits(c::NamedTuple)
    return (
        hit_x = [c.hit_x[2], c.hit_x[1]], 
        hit_y = [c.hit_y[2], c.hit_y[1]], 
        hit_z = [c.hit_z[2], c.hit_z[1]], 
        hit_edep = [c.hit_edep[2], c.hit_edep[1]]
        )
end

@inline function is_valid_2hit(evt)
    length(evt.hit_x) == 2 && hypot(evt.hit_x[2] - evt.hit_x[1],
        evt.hit_y[2] - evt.hit_y[1],
        evt.hit_z[2] - evt.hit_z[1]) > 3.0u"mm"
end

function getz(file; name="segBEGe", center=81.76361317572471, ew = 8.0u"keV")
    icpc, czt = LHDataStore(file) do lhd
        lhd[name][:], lhd["czt"][:]
    end
    hv = getV(file)
    ec = econv[hv]
    icpc_e = ec*icpc.DAQ_energy * u"keV"
    czt_e = uconvert.(u"keV", (sum.(czt.hit_edep)))
    idx = intersect(findall(x -> abs(x - Cs_energy) ≤ ew, icpc_e+czt_e), 
                    findall(x -> 250.0u"keV" ≤ x ≤ 440.0u"keV", icpc_e))
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


"""
    get_all_z(sourcedir; center=CNTR_VALUE)

loop through all files in `sourcedir` and return for each file the radius, 
measuretime and reconstructed z's
"""
function get_all_z(sourcedir; name="segBEGe", center=81.76361317572471)
    ffiles = filter(x -> endswith(x, "preprocessed.lh5"), readdir(sourcedir))
    R = zeros(length(ffiles))
    mtime = zeros(length(ffiles))
    z = Vector{Float64}[]
    for i=eachindex(ffiles)
        R[i] = getR(ffiles[i])
        mtime[i] = getM(ffiles[i])
        rec_zs = getz(joinpath(sourcedir, ffiles[i]); name, center, ew = 20.0u"keV")
        push!(z, rec_zs[:, 1])
    end
    mtime, R, z
end

function get_z_and_wavefroms(file, hv; i=1, name="segBEGe", ew= 8.0u"keV", Δz=1)
    det, czt = LHDataStore(file) do lhd
        lhd[name][:], lhd["czt"][:]
    end
    det = begin
        idx = findall(x -> x == i, det.chid)
        det[idx]
    end
    println("here")
    ec = econv[hv]
    icpc_e = ec*det.DAQ_energy * u"keV"
    czt_e = uconvert.(u"keV", (sum.(czt.hit_edep)))
    idx = intersect(findall(x -> abs(x - Cs_energy) ≤ ew, icpc_e+czt_e), 
                    findall(x -> 250.0u"keV" ≤ x ≤ 440.0u"keV", icpc_e))
    det_hits = view(det, idx)
    czt_hits = view(czt, idx)
    R = cntr - getR(file)

    @info "Reconstructing Z from two hit events at R = $R"
    idx_2h = findall(is_valid_2hit, czt_hits);
    idx_1h = findall(x -> x == 1, czt_hits.evt_nhits)
    det_1hit = view(det_hits, idx_1h)
    czt_1hit = view(czt_hits, idx_1h)
    det_2hit = view(det_hits, idx_2h)
    czt_2hit = view(czt_hits, idx_2h)

    @time z_rec_1hit = get_z_from_energies.(det_1hit, czt_1hit, R, hv)
    @time z_rec_2hit = get_z_from_2_hit_events.(det_2hit, czt_2hit, R; Δz, hv)
    idx_val = findall(x -> x[1] != 0, z_rec_2hit)
    z_rec_2hit_val = [x[2] for x in z_rec_2hit[idx_val]]
    return z_rec_2hit_val, z_rec_1hit, det_2hit.samples[idx_val], det_hits.samples[idx_1h] 
end

export get_z_and_wavefroms

function normalizewf(x::Vector{T}, l::Int) where {T}
    max = mean(x[end-l:end])
    x/max
end

function normalizewf(x::Vector{Vector{T}}, l::Int) where {T}
    max = sum([mean(x[i][end-l:end]) for i=eachindex(x)])
    sum([wfm/max for wfm=wfs_aligned])
end

function reconstruct_at_radius(file, hv; Δz=1, window=(500, 500), τ=51.8, 
χ2_max=3, l1=300)
    z_rec_2h, z_rec_1h, wf_2hit, wf_1hit = get_z_and_wavefroms(file, hv)
    superpulses_Cs = []
    z_Cs = collect(0:Δz:40)
    mask = zeros(Bool, length(z_Cs))
    for i=eachindex(z_Cs)
        # filter waveforms from validated two hit events
        idxz2 = findall(z -> abs(z - z_Cs[i]) < Δz, z_rec_2h)
        length(wfs_at_z) == 0 && break
        wfs2_at_z = baseline_corr.(wf_2hit[idxz2])
        # TODO check usefulleness of is_singlesite
        # wfs_at_z = filter(wf -> is_singlesite(wf), wfs_at_z)
        # sprpls = superpulse(wfs_at_z)
        wlength = sum(window) + 1
        wfs2_at_z = decay_correction.(wfs2_at_z, exp(-0.004/τ))
        wfs_aligned = time_align.(wfs2_at_z; 0.5, window)
        filter!(wfm -> length(wfm) == wlength, wfs_aligned)
        length(wfs_aligned) == 0 && break

        # wmax = sum([mean(wfm[end-l1:end]) for wfm=wfs_aligned])
        # sp = sum([wfm/wmax for wfm=wfs_aligned])
        sp = normalizewf(wfs_aligned, l1)
        χ2 = [chi_sq_wfs(normalizewf(wfm, l1), sp) for wfm=wfs_aligned]
        idx_chi2 = findall(x -> x < χ2_max, χ2)
        length(idx_chi2) == 0 && break

        # wmax = sum([mean(wfm[end-l1:end]) for wfm=wfs_aligned[idx_chi2]])
        # sp = sum([wfm/wmax for wfm=wfs_aligned])
        sp = normalizewf(wfs_aligned[idx_chi2], l1)

        idxz1 = findall(z -> abs(z - z_Cs[i]) < Δz, z_rec_1h)
        wfs1_at_z = baseline_corr.(wf_1hit[idxz1])
        wfs1_at_z = decay_correction.(wfs1_at_z, exp(-0.004/τ))
        wfs1_aligned = time_align.(wfs1_at_z; 0.5, window)

    end
end

export reconstruct_at_radius

baseline_corr(wf; m=500) = wf .- mean(wf[1:m])

# function superpulse(wfs::Vector{T}; τ::Float64=51.8, 
# window::Tuple{Int, Int}=(500, 500)) where {T}
#     wlength = sum(window) + 1
#     wfs = decay_correction.(wfs, exp(-0.004/τ))
#     wfs = time_align.(wfs, at=0.5, window=(500, 500))
#     filter!(wf -> length(wf) == wlength, wfs)
#     if length(wfs) > 0
#         w = sum([mean(wf[end-300:end]) for wf in wfs])
#         return sum([wf/w for wf in wfs])
#     else
#         return fill(NaN, sum(window) + 1)
#     end
# end

# export superpulse


"""
    time_align(wf::AbstractVector{T}; at::T, window::Tuple{Int, Int})
"""
@fastmath function time_align(wf::AbstractVector{T}; at::T=0.5, 
window::Tuple{Int, Int}=(500, 500))::Vector{T} where {T<:AbstractFloat}
    @inbounds begin
        wfmax = mean(wf[end-300:end])
        idx = Int(find_crossing(smooth20(wf), at*wfmax, interpolate = false))
        if idx ≤ window[1] || length(wf) - window[2] < idx
            @warn "Time outside of bounds, wf not aligned"
            return wf
        end
        wf[idx-window[1]:idx+window[2]]
    end
end

@fastmath function decay_correction(wf::AbstractVector{T}, decay_factor::T
)::Vector{T} where {T<:AbstractFloat}
    rp = Vector{T}(undef, length(wf))
    rp[1] = wf[1]
    @inbounds for i in 2:length(rp)
        rp[i] = wf[i] + rp[i-1] - wf[i-1] * decay_factor
    end
    rp
end

function smooth20(wv::Vector{T})::Vector{T} where {T <: AbstractFloat}
    w::Vector{T} = T[ 0.0,
                     0.010329700377119983,
                     0.022145353094771694,
                     0.034827599769728275,
                     0.04766125969628491,
                     0.05988618988841385,
                     0.07075302574183234,
                     0.07957926681837389,
                     0.08580108206476868,
                     0.0890165225487064,
                     0.0890165225487064,
                     0.08580108206476868,
                     0.07957926681837389,
                     0.07075302574183234,
                     0.05988618988841385,
                     0.04766125969628491,
                     0.034827599769728275,
                     0.022145353094771694,
                     0.010329700377119983,
                     0.0 ]
    DSP.filtfilt(w, T[1], wv)
end

@fastmath function find_crossing(wf::AbstractArray{T, 1}, thold::T; window::Tuple{Int64, Int64} = (100, 1100), interpolate::Bool = true)::T where {T<:AbstractFloat}
    @inbounds begin
        cross = -1
        i = window[1]
        while cross == -1 && i < window[2]
            if wf[i] ≤ thold && wf[i+1] ≥ thold
                cross = i
            else
                i = i+1
            end
        end
        if cross != -1 && interpolate
            cross = cross + (thold - wf[cross])/(wf[cross+1] - wf[cross])     
        end
        return cross
    end
end

function chi_sq_wfs(x, y, window=400:600, l=150)
    varx = var(x[1:l])
    vary = var(y[1:l])
    sum((x[window] - y[window]).^2) / ((length(window) - 1)*(varx + vary))
end