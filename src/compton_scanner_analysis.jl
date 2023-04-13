# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

function reconstruct_z(file::AbstractString; name::AbstractString = "segBEGe", 
center::QuantityMM{Float64} = cntr, ΔE = 8.0u"keV", Δz::QuantityMM = 2.0u"mm", 
I_core::ClosedInterval{<:RealQuantity} = 250.0u"keV"..440.0u"keV")
    
    det::detTable, czt::cztTable, idx_c::Int, ec::typeof(Cs_energy) = 
        read_preprocessed_file(file, name)
    core = @view det[findall(det.chid .== idx_c)]
    hv::Float64 = ustrip(getV(file))
    core_e = ec*core.DAQ_energy
    E_tot = uconvert.(u"keV", (sum.(czt.hit_edep))) + core_e
    I_tot = Cs_energy ± ΔE
    idx = intersect(findall(in(I_tot), E_tot), findall(in(I_core), core_e))
    core_hits = view(core, idx)
    czt_hits = view(czt, idx)
    R::QuantityMM{Float64} = center - getR(file)

    @info "Reconstructing Z from two hit events at R = $R"
    idx_2h = findall(is_valid_2hit, czt_hits)
    czt_2hit = view(czt_hits, idx_2h)
    det_2hit = view(core_hits, idx_2h)
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
    get_all_z(sourcedir::AbstractString; name::AbstractString="segBEGe", 
    center::QuantityMM{Float64}=CNTR_VALUE, ΔE = 8.0u"keV")

loop through all files in `sourcedir` and return for each file the radius, 
measuretime and reconstructed z's
"""
function get_all_z(sourcedir::AbstractString; name::AbstractString="segBEGe", 
center::QuantityMM{Float64} = cntr, ΔE = 8.0u"keV")

    ffiles = filter(x -> endswith(x, "preprocessed.lh5"), readdir(sourcedir))
    R = zeros(QuantityMM{Float64}, length(ffiles))
    mtime = zeros(typeof(1.0u"s"), length(ffiles))
    z = Vector{Vector{QuantityMM{Float64}}}(undef, length(ffiles))
    for i=eachindex(ffiles)
        file = joinpath(sourcedir, ffiles[i])
        R[i] = getR(file)
        mtime[i] = getM(file)
        rec_zs = reconstruct_z(file; name, center, ΔE = ΔE)
        z[i] = rec_zs[:, 1]
    end
    ustrip.(mtime), ustrip.(R), broadcast(z -> ustrip.(z), z) # remove units for now
end

"""
    get_z_and_waveforms(file::AbstractString, hv::Float64, chid::Int;
    center::QuantityMM{Float64} = cntr, name::AbstractString = "segBEGe",
    ΔE = 8.0u"keV", Δz::QuantityMM = 2.0u"mm", 
    I_core::ClosedInterval{<:RealQuantity} = 250.0u"keV"..440.0u"keV")
    
For a specified preprocessed `file`, read in the core and the segment
table (given by `chid`), compute the reconstructed z values of one and 
two hit events and return those together with the corresponding waveforms
of the specified segment.
"""
function get_z_and_waveforms(file::AbstractString, hv::Float64, chid::Int; 
center::QuantityMM{Float64} = cntr, name::AbstractString = "segBEGe", 
ΔE = 8.0u"keV", Δz::QuantityMM = 2.0u"mm", 
I_core::ClosedInterval{<:RealQuantity} = 250.0u"keV"..440.0u"keV")

    det::detTable, czt::cztTable, idx_c::Int, ec::typeof(Cs_energy) = 
        read_preprocessed_file(file, name)
    core = view(det, findall(det.chid .== idx_c))
    seg = view(det, findall(det.chid .== chid))
    @assert length(core) == length(czt) "core and czt do not have the same number of events"
    @assert length(seg) == length(czt) "segment $chid and czt do not have the same number of events"
    core_e = ec * core.DAQ_energy
    E_tot = uconvert.(u"keV", (sum.(czt.hit_edep))) + core_e
    I_tot = Cs_energy ± ΔE
    idx = intersect(findall(in(I_tot), E_tot), findall(in(I_core), core_e))
    core_hits = view(core, idx)
    seg_hits = view(seg, idx)
    czt_hits = view(czt, idx)
    R::QuantityMM{Float64} = center - getR(file)

    idx_2h = findall(is_valid_2hit, czt_hits)
    core_2hit = view(core_hits, idx_2h)
    seg_2hit = view(seg_hits, idx_2h)
    czt_2hit = view(czt_hits, idx_2h)
    z_rec_2hit = get_z_from_2_hit_events.(core_2hit, czt_2hit, R; Δz, hv)
    idx_val = findall(x -> x[1] != 0, z_rec_2hit)
    z_rec_2hit_val = [x[2] for x in z_rec_2hit[idx_val]]

    idx_1h = findall(czt_hits.evt_nhits .== 1)
    core_1hit = view(core_hits, idx_1h)
    seg_1hit = view(seg_hits, idx_1h)
    czt_1hit = view(czt_hits, idx_1h)
    z_rec_1hit = get_z_from_energies.(core_1hit, czt_1hit, R, hv)

    return z_rec_2hit_val, z_rec_1hit, view(seg_2hit, idx_val).samples, 
        seg_1hit.samples
end


# TODO: find better solution for handling kwargs
"""
    reconstruct_at_radius(file::AbstractString), hv::Float64, chid::Int; 
    name::AbstractString = "segBEGe", Δz::TT=2.0u"mm", zbin::TT = 2.0u"mm", 
    window::Tuple{Int, Int}=(500, 500), τ::Float64=51.8, χ2_max::Float64=3., 
    l1::Int=300, ΔE = 8.0u"keV", baseline_length::Int = 500, 
    verbose::Bool = true) where {TT <: QuantityMM}

Given a `file`, a specified voltage `hv` and the corresponding channel
id `chid`, build the superpulses from two hit and one hit validated hits.  
"""
function reconstruct_at_radius(file::AbstractString, hv::Float64, chid::Int; 
name::AbstractString = "segBEGe", Δz::TT = 2.0u"mm", zbin::TT = 2.0u"mm",
window::Tuple{Int, Int}=(500, 500), τ::Float64=51.8, χ2_max::Float64 = 3.,
l1::Int = 300, ΔE = 8.0u"keV", baseline_length::Int = 500, 
verbose::Bool = true) where {TT <: QuantityMM}    

    z2h, z1h, wf2h, wf1h = get_z_and_waveforms(file, hv, chid, name = name, ΔE = ΔE, Δz = Δz)
    L = sum(window)+1
    z_Cs = collect(zero(TT):zbin:TT(40u"mm"))
    superpulses_Cs = nestedview(Matrix{Float64}(undef, L, length(z_Cs)))
    mask = zeros(Bool, length(z_Cs))
    if verbose 
        prog = ProgressUnknown("Reconstructing superpulses for channel $chid \
        at z =") 
    end
    for i=eachindex(z_Cs)
        verbose && ProgressMeter.update!(prog, i)
        I_z1 = z_Cs[i] ± zbin
        I_z2 = z_Cs[i] ± zbin/2

        # select all two-hit events in the given z-bin
        idxz2 = findall(in(I_z2), z2h)
        length(idxz2) == 0 && continue

        # correct the pulses for baseline and τ, and time-align them
        # TODO check usefulness of is_singlesite for two-hit events
        wf2h_at_z = baseline_corr!.(wf2h[idxz2]; m=baseline_length)
        wf2h_at_z = decay_correction.(wf2h_at_z, exp(-0.004/τ))
        wfs_aligned = time_align.(wf2h_at_z; p=0.5, window=window, l=l1)
        filter!(wfm -> length(wfm) == L, wfs_aligned)
        length(wfs_aligned) == 0 && continue

        # determine a two-hit superpulse and discard obvious outliers
        sp = normalizewf(wfs_aligned; tail_length=l1)
        χ2 = [chi_sq_wfs(normalizewf!(wfm; tail_length=l1), sp) for wfm=wfs_aligned]
        chi_2h = findall(χ2 .< χ2_max)
        length(chi_2h) == 0 && continue

        # create superpulse from two-hit waveforms which passed the first χ2 cut
        sp = normalizewf(wfs_aligned[chi_2h]; tail_length=l1)

        # select all one-hit events in the given z-bin
        idxz1 = findall(in(I_z1), z1h)

        # correct the pulses for baseline and τ, normalize and time-align them
        wfs1_at_z = baseline_corr!.(wf1h[idxz1]; m=baseline_length)
        wfs1_at_z = decay_correction.(wfs1_at_z, exp(-0.004/τ))
        wfs1_aligned = time_align.(wfs1_at_z; p=0.5, window=window, l=l1)
        filter!(wfm -> length(wfm) == L, wfs1_aligned)
        if length(wfs1_aligned) == 0
            superpulses_Cs[i] .= sp
            mask[i] = true
            continue
        end

        # select all one-hit events with χ2 < χ2_max compared to the two-hit superpulse
        χ1h = [chi_sq_wfs(normalizewf!(wfm; tail_length=l1), sp) for wfm=wfs1_aligned]
        chi_1h = findall(χ1h .< χ2_max)
        sp = normalizewf(
            vcat(wfs_aligned[chi_2h], wfs1_aligned[chi_1h]); tail_length=l1)
        superpulses_Cs[i] .= sp
        mask[i] = true
    end
    superpulses_Cs, mask, ustrip.(z_Cs)
end