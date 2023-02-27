# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

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
    idx_c::Int = 1, Δz::TT = 2.0u"mm", zbin::TT = 2.0u"mm",
    window::Tuple{Int, Int}=(500, 500), τ::Float64=51.8, χ2_max::Float64=3.,
    l1::Int=300, ew = 8.0u"keV", baseline_samples::Int = 500, verbose::Bool = true) where {TT <: QuantityMM}
    
    z2h, z1h, wf2h, wf1h = get_z_and_waveforms(file, hv, name = name, idx_c = idx_c, ew = ew, Δz = Δz)
    wlength = sum(window)+1
    z_Cs = collect(zero(TT):zbin:TT(40u"mm"))
    superpulses_Cs = [zeros(Float64, wlength) for _ in eachindex(z_Cs)]
    mask = zeros(Bool, length(z_Cs))
    if verbose prog = ProgressUnknown("Reconstructing superpulses at z =") end
    for i=eachindex(z_Cs)
        verbose && ProgressMeter.update!(prog, i)

        # select all two-hit events in the given z-bin
        idxz2 = findall(z -> abs(z - z_Cs[i]) < zbin/2, z2h)
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