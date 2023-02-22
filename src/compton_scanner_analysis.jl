in_mm(x::Number) = ustrip(uconvert(u"mm", x))

function get_z_from_2_hit_events(
        det::detTuple, czt::cztTuple, R::QuantityMM{Float64}; 
        Δz::QuantityMM = 2.0u"mm", hv::Float64
    )::Tuple{Int, QuantityMM{Float64}, QuantityMM{Float64}}

    #if econv*det.DAQ_energy < 250
        #return (false, -1)
    #end
    #get z coordinates for scatter point in segBEGe
    zθ::QuantityMM{Float64} = get_z_from_energies(det, czt, R, hv)
    zα::Vector{QuantityMM{Float64}} = get_z_from_camera(czt, R)

    #compare them: if they agree within Δz, keep them
    for z in zα
        if abs(z - zθ) < Δz 
            return (1, zθ, z)
        end
    end
    
    #swaphits and try again
    czts = swap_CZT_hits(czt)
    zθ = get_z_from_energies(det, czts, R, hv)
    zα = get_z_from_camera(czts, R)
    for z in zα
        if abs(z - zθ) < Δz
            return (2, zθ, z)
        end
    end

    #else discard them
    return (0, NaN*u"mm", NaN*u"mm")
end


function get_z_from_energies(det::detTuple, czt::cztTuple, R::QuantityMM{Float64}, hv::Float64)::QuantityMM{Float64}
    T = QuantityMM{Float64}
    # TODO: add flag if we want to rely more on camera or DAQ_energies?
    # if ge not depleted -> worse energy resolution -> rely more on camera 
    # energy resolution
    # DAQ_energies
    # cf = econv[hv]
    # θ = compton_angle(cf*det.DAQ_energy*u"keV"+sum(czt.hit_edep), sum(czt.hit_edep))
    # CZT energies
    θ::Float64 = compton_angle(Cs_energy, sum(czt.hit_edep)) 
    T(czt.hit_z[1]) + hypot(T(czt.hit_x[1]), T(czt.hit_y[1] - R)) * cot(θ)
end


function get_z_from_camera(czt::cztTuple, R::QuantityMM{Float64})
    let x = czt.hit_x, y = czt.hit_y, z = czt.hit_z
        α::Float64 = compton_angle(sum(czt.hit_edep), czt.hit_edep[2])
        zα = []

        if !(isnan(α))
            try
                #define the cone and get intersections z1,z2 with the beam axis
                camhit1::Vector{QuantityMM{Float64}} = [x[1], y[1], z[1]]
                camhit2::Vector{QuantityMM{Float64}} = [x[2], y[2], z[2]]
                cone = Cone(camhit1, camhit1-camhit2, α)
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
end


function get_possible_z_from_camera(cone::Cone{T,TT}, R::TT) where {T <: AbstractFloat, TT <: QuantityMM{T}}

    #get relevant cone parameters
    u = unit(TT)
    H = cone.axis
    c0 = cone.origin - [zero(TT), R, zero(TT)]
    cosα = cos(cone.α)

    #solve quadartic equation for z
    A = H[3]^2 - cosα^2
    B = -2*H[3]*dot(H,c0) + 2*c0[3]*cosα^2
    C = dot(H,c0)^2 - norm(c0)^2 * cosα^2
    z1 = (-B - sqrt(B^2 - 4*A*C)) / (2*A)
    z2 = (-B + sqrt(B^2 - 4*A*C)) / (2*A)

    return z1, z2

end


function validate_z(z::TT, cone::Cone{T,TT}, R::TT; Δα::T = T(1e-6)) where {T <: AbstractFloat, TT <: QuantityMM{T}}

    #calculate the angle from the vectors
    tmp = TT[zero(TT), R, z] - cone.origin
    αnew = acos(dot(cone.axis, tmp/norm(tmp)))

    #keep only true angles
    return abs(cone.α - αnew) < Δα
end

function swap_CZT_hits(czt_old::cztTuple)::cztTuple
    # TODO: check if the deepcopy is needed or whether
    # we can reverse the Vectors in the original cztTuple
    czt = deepcopy(czt_old)
    reverse!(czt.hit_x)
    reverse!(czt.hit_y)
    reverse!(czt.hit_z)
    reverse!(czt.hit_edep)
    return czt
end

@inline function is_valid_2hit(czt::cztTuple)::Bool
    length(czt.hit_x) == 2 && hypot(czt.hit_x[2] - czt.hit_x[1],
        czt.hit_y[2] - czt.hit_y[1],
        czt.hit_z[2] - czt.hit_z[1]) > 3.0u"mm"
end

function reconstruct_z(file::AbstractString; name::AbstractString = "segBEGe", 
    center::QuantityMM{Float64} = cntr, ew = 8.0u"keV", Δz::QuantityMM = 2.0u"mm")
    det::detTable, czt::cztTable = read_preprocessed_file(file, name)
    hv::Float64 = ustrip(getV(file))
    ec = econv[hv] * u"keV"
    det_e = ec*det.DAQ_energy
    czt_e = uconvert.(u"keV", (sum.(czt.hit_edep)))
    idx = intersect(findall(x -> abs(x - Cs_energy) ≤ ew, det_e+czt_e), 
                    findall(x -> 250.0u"keV" ≤ x ≤ 440.0u"keV", det_e))
    det_hits = view(det, idx)
    czt_hits = view(czt, idx)
    R::QuantityMM{Float64} = center - getR(file)

    @info "Reconstructing Z from two hit events at R = $R"
    idx_2h = findall(is_valid_2hit, czt_hits)
    czt_2hit = view(czt_hits, idx_2h)
    det_2hit = view(det_hits, idx_2h)
    # TODO: resolve allocation issues by passing fixed empty array
    zrec2hit = get_z_from_2_hit_events.(det_2hit, czt_2hit, R; Δz, hv)
    idx_val_1 = findall(x -> x[1] == 1, zrec2hit)
    idx_val_2 = findall(x -> x[1] == 2, zrec2hit)
    idx_val = vcat(idx_val_1, idx_val_2)

    # TODO: decide on what reconstructed z we really want 
    # from energies (θ)    -> 2
    # from camera cone (α) -> 3
    vcat([[x[2] x[3]] for x in view(zrec2hit, idx_val)]...)
end


"""
    get_all_z(sourcedir; center=CNTR_VALUE)

loop through all files in `sourcedir` and return for each file the radius, 
measuretime and reconstructed z's
"""
function get_all_z(sourcedir::AbstractString; name::AbstractString="segBEGe", 
    center::QuantityMM{Float64} = cntr, ew = 8.0u"keV")
    ffiles = filter(x -> endswith(x, "preprocessed.lh5"), readdir(sourcedir))
    R = zeros(QuantityMM{Float64}, length(ffiles))
    mtime = zeros(typeof(1.0u"s"), length(ffiles))
    z = Vector{QuantityMM{Float64}}[]
    for i=eachindex(ffiles)
        R[i] = getR(ffiles[i])
        mtime[i] = getM(ffiles[i])
        rec_zs = reconstruct_z(joinpath(sourcedir, ffiles[i]); name, center, ew = ew)
        push!(z, rec_zs[:, 1])
    end
    ustrip.(mtime), ustrip.(R), broadcast(z -> ustrip.(z), z) # remove units for now
end

function get_z_and_waveforms(
        file::AbstractString, hv::Float64; center::QuantityMM{Float64} = cntr, idx_c::Int = 1, 
        name::AbstractString = "segBEGe", ew = 8.0u"keV", Δz::QuantityMM = 2.0u"mm"
    )

    det, czt = read_preprocessed_file(file, name)
    det = det[findall(x -> x == idx_c, det.chid)]
    @assert length(det) == length(czt) "$name and czt do not have the same number of events"
    ec = econv[hv] * u"keV"
    det_e = ec*det.DAQ_energy
    czt_e = uconvert.(u"keV", (sum.(czt.hit_edep)))
    idx = intersect(findall(x -> abs(x - Cs_energy) ≤ ew, det_e+czt_e), 
                    findall(x -> 250.0u"keV" ≤ x ≤ 440.0u"keV", det_e))
    det_hits = view(det, idx)
    czt_hits = view(czt, idx)
    R::QuantityMM{Float64} = center - getR(file)

    # @info "Reconstructing Z from two hit events at R = $R"
    idx_2h = findall(is_valid_2hit, czt_hits)
    idx_1h = findall(x -> x == 1, czt_hits.evt_nhits)
    det_1hit = view(det_hits, idx_1h)
    czt_1hit = view(czt_hits, idx_1h)
    det_2hit = view(det_hits, idx_2h)
    czt_2hit = view(czt_hits, idx_2h)

    z_rec_1hit = get_z_from_energies.(det_1hit, czt_1hit, R, hv)
    z_rec_2hit = get_z_from_2_hit_events.(det_2hit, czt_2hit, R; Δz, hv)
    idx_val = findall(x -> x[1] != 0, z_rec_2hit)
    z_rec_2hit_val = [x[2] for x in z_rec_2hit[idx_val]]
    return z_rec_2hit_val, z_rec_1hit, det_2hit.samples[idx_val], det_hits.samples[idx_1h] 
end


# TODO: find better solution for handling kwargs
"""
    reconstruct_at_radius(file::AbstractString), hv::Float64; Δz::QuantityMM=2.0u"mm",
    window::Tuple{Int, Int}=(500, 500), τ::Float64=51.8,
    χ2_max::Float64=3., l1::Int=300, ew = 8.0u"keV")

Given a `file` and a specified voltage `hv`, build the superpulses from 
two hit and one hit validated hits.  
"""
function reconstruct_at_radius(file::AbstractString, hv::Float64; name::AbstractString = "segBEGe", 
    idx_c::Int = 1, Δz::TT = 2.0u"mm", zbin::TT = 1.0u"mm",
    window::Tuple{Int, Int}=(500, 500), τ::Float64=51.8, χ2_max::Float64=3.,
    l1::Int=300, ew = 8.0u"keV", baseline_samples::Int = 500, verbose::Bool = true) where {TT <: QuantityMM}
    
    z2h, z1h, wf2h, wf1h = get_z_and_waveforms(file, hv, name = name, idx_c = idx_c, ew = ew, Δz = Δz)
    wlength = sum(window)+1
    z_Cs = collect(zero(TT):2*zbin:TT(40u"mm"))
    superpulses_Cs = [zeros(Float64, wlength) for _ in eachindex(z_Cs)]
    mask = zeros(Bool, length(z_Cs))
    if verbose prog = ProgressUnknown("Reconstructing superpulses at z =") end
    for i=eachindex(z_Cs)
        verbose && ProgressMeter.update!(prog, i)

        # select all two-hit events in the given z-bin
        idxz2 = findall(z -> abs(z - z_Cs[i]) < zbin, z2h)
        length(idxz2) == 0 && continue

        # correct the pulses for baseline and τ, and time-align them
        # TODO check usefulness of is_singlesite for two-hit events
        wf2h_at_z = baseline_corr.(wf2h[idxz2]; m=baseline_samples)
        wf2h_at_z = decay_correction.(wf2h_at_z, exp(-0.004/τ))
        wfs_aligned = time_align.(wf2h_at_z; p=0.5, window=window, l=l1)
        filter!(wfm -> length(wfm) == wlength, wfs_aligned)
        length(wfs_aligned) == 0 && continue

        # determine a two-hit superpulse and discard obvious outliers
        sp = normalizewf(wfs_aligned; l=l1)
        χ2 = [chi_sq_wfs(normalizewf(wfm; l=l1), sp) for wfm=wfs_aligned]
        chi_2h = findall(x -> x < χ2_max, χ2)
        length(chi_2h) == 0 && continue

        # create superpulse from two-hit waveforms which passed the first χ2 cut
        sp = normalizewf(wfs_aligned[chi_2h]; l=l1)

        # select all one-hit events in the given z-bin
        idxz1 = findall(z -> abs(z - z_Cs[i]) < zbin, z1h)

        # correct the pulses for baseline and τ, normalize and time-align them
        wfs1_at_z = baseline_corr.(wf1h[idxz1]; m=baseline_samples)
        wfs1_at_z = decay_correction.(wfs1_at_z, exp(-0.004/τ))
        wfs1_aligned = time_align.(wfs1_at_z; p=0.5, window=window, l=l1)
        filter!(wfm -> length(wfm) == wlength, wfs1_aligned)
        if length(wfs1_aligned) == 0
            superpulses_Cs[i] = sp
            mask[i] = true
            continue
        end

        # select all one-hit events with χ2 < χ2_max compared to the two-hit superpulse
        χ1h = [chi_sq_wfs(normalizewf(wfm; l=l1), sp) for wfm=wfs1_aligned]
        chi_1h = findall(χ -> χ < χ2_max, χ1h)
        sp = normalizewf(vcat(wfs_aligned[chi_2h], wfs1_aligned[chi_1h]); l=l1)
        superpulses_Cs[i] = sp
        mask[i] = true
    end
    superpulses_Cs, mask, ustrip.(z_Cs)
end