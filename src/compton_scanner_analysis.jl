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
    T = typeof(1.0*unit(eltype(c.hit_x)))
    # TODO: add flag if we want to rely more on camera or DAQ_energies?
    # if ge not depleted -> worse energy resolution -> rely more on camera 
    # energy resolution
    # DAQ_energies
    # cf = econv[hv]
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
    return abs(cone.α - αnew) < Δα
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

function getz(file; name="segBEGe", center=cntr, ew = 8.0u"keV", Δz=2)
    det, czt = LHDataStore(file) do lhd
        lhd[name][:], lhd["czt"][:]
    end
    hv = getV(file)
    ec = econv[hv] * u"keV"
    det_e = ec*det.DAQ_energy
    czt_e = uconvert.(u"keV", (sum.(czt.hit_edep)))
    idx = intersect(findall(x -> abs(x - Cs_energy) ≤ ew, det_e+czt_e), 
                    findall(x -> 250.0u"keV" ≤ x ≤ 440.0u"keV", det_e))
    det_hits = view(det, idx)
    czt_hits = view(czt, idx)
    R = center - getR(file)

    @info "Reconstructing Z from two hit events at R = $R"
    idx_2h = findall(is_valid_2hit, czt_hits)
    czt_2hit = view(czt_hits, idx_2h)
    det_2hit = view(det_hits, idx_2h)
    # TODO: resolve allocation issues by passing fixed empty array
    zrec2hit = get_z_from_2_hit_events.(det_2hit, czt_2hit, R; Δz = Δz, hv)
    idx_val_1 = findall(x -> x[1] == 1, zrec2hit)
    idx_val_2 = findall(x -> x[1] == 2, zrec2hit)
    idx_val = vcat(idx_val_1, idx_val_2)

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
function get_all_z(sourcedir::AbstractString; name::AbstractString="segBEGe", center=cntr, ew = 8.0u"keV")
    ffiles = filter(x -> endswith(x, "preprocessed.lh5"), readdir(sourcedir))
    R = zeros(length(ffiles))
    mtime = zeros(length(ffiles))
    z = Vector{Float64}[]
    for i=eachindex(ffiles)
        R[i] = getR(ffiles[i])
        mtime[i] = getM(ffiles[i])
        rec_zs = getz(joinpath(sourcedir, ffiles[i]); name, center, ew = ew)
        push!(z, rec_zs[:, 1])
    end
    mtime, R, z
end

function get_z_and_waveforms(file::AbstractString, hv; center = cntr, idx_c=1, name="segBEGe", ew= 8.0u"keV", Δz=1)
    det, czt = LHDataStore(file) do lhd
        lhd[name][:], lhd["czt"][:]
    end
    det = det[findall(x -> x == idx_c, det.chid)]
    @assert length(det) == length(czt) "$name and czt do not have the same number of events"
    ec = econv[hv] * u"keV"
    det_e = ec*det.DAQ_energy
    czt_e = uconvert.(u"keV", (sum.(czt.hit_edep)))
    idx = intersect(findall(x -> abs(x - Cs_energy) ≤ ew, det_e+czt_e), 
                    findall(x -> 250.0u"keV" ≤ x ≤ 440.0u"keV", det_e))
    det_hits = view(det, idx)
    czt_hits = view(czt, idx)
    R = center - getR(file)

    # @info "Reconstructing Z from two hit events at R = $R"
    idx_2h = findall(is_valid_2hit, czt_hits);
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
    reconstruct_at_radius(file::AbstractString), hv::Float64; Δz::Int=1,
    window::Tuple{Int, Int}=(500, 500), τ::Float64=51.8,
    χ2_max::Float64=3., l1::Int=300, ew = 8.0u"keV")

Given a `file` and a specified voltage `hv`, build the superpulses from 
two hit and one hit validated hits.  
"""
function reconstruct_at_radius(file::AbstractString, hv::Float64; Δz::Int=1, 
window::Tuple{Int, Int}=(500, 500), τ::Float64=51.8, χ2_max::Float64=3., 
l1::Int=300, ew = 8.0u"keV", baseline_samples::Int = 500, verbose::Bool = true)
    z2h, z1h, wf2h, wf1h = get_z_and_waveforms(file, hv, ew = ew)
    wlength = sum(window)+1
    superpulses_Cs = []
    z_Cs = collect(0:2:40)
    mask = zeros(Bool, length(z_Cs))
    if verbose prog = ProgressUnknown("Reconstructing superpulses at z =") end
    for i=eachindex(z_Cs)
        verbose && ProgressMeter.update!(p, i)
        idxz2 = findall(z -> abs(z - z_Cs[i]) < Δz, z2h)
        length(idxz2) == 0 && break
        # TODO check usefulness of is_singlesite
        wf2h_at_z = baseline_corr.(wf2h[idxz2]; m=baseline_samples)
        wf2h_at_z = decay_correction.(wf2h_at_z, exp(-0.004/τ))
        wfs_aligned = time_align.(wf2h_at_z; p=0.5, window=window, l=l1)
        filter!(wfm -> length(wfm) == wlength, wfs_aligned)
        length(wfs_aligned) == 0 && break

        sp = normalizewf(wfs_aligned; l=l1)
        χ2 = [chi_sq_wfs(normalizewf(wfm; l=l1), sp) for wfm=wfs_aligned]
        chi_2h = findall(x -> x < χ2_max, χ2)
        length(chi_2h) == 0 && break

        # new superpulse from 2hit waveforms which passed first χ2
        sp = normalizewf(wfs_aligned[chi_2h]; l=l1)
        idxz1 = findall(z -> abs(z - z_Cs[i]) < Δz, z1h)
        wfs1_at_z = baseline_corr.(wf1h[idxz1]; m=baseline_samples)
        wfs1_at_z = decay_correction.(wfs1_at_z, exp(-0.004/τ))
        wfs1_aligned = time_align.(wfs1_at_z; p=0.5, window=window, l=l1)
        filter!(wfm -> length(wfm) == wlength, wfs1_aligned)
        if length(wfs1_aligned) == 0
            append!(superpulses_Cs, [sp])
            mask[i] = true
            break
        end

        χ1h = [chi_sq_wfs(normalizewf(wfm; l=l1), sp) for wfm=wfs1_aligned]
        chi_1h = findall(χ -> χ < χ2_max, χ1h)
        sp = normalizewf(vcat(wfs_aligned[chi_2h], wfs1_aligned[chi_1h]); l=l1)
        append!(superpulses_Cs, [sp])
        mask[i] = true
    end
    superpulses_Cs, mask, z_Cs
end

function normalizewf(x; l::Int=300)
    max = mean(x[end-l:end])
    x/max
end

function normalizewf(x::AbstractVector{<:AbstractVector{T}}; l::Int=300
) where {T}
    max = sum([mean(x[i][end-l:end]) for i=eachindex(x)])
    sum([x[i]/max for i=eachindex(x)])
end

# TODO: check if ! functions would give better performance 
# (certainly better allocs)
# baseline_corr!(x::AbstractVector{T}; m::Int=500) where {T} = begin
#     x̅ = mean(x[1:m])
#     x .-= x̅
# end

baseline_corr(x::AbstractVector{T}; m::Int=500) where {T} = begin
    x̅ = mean(x[1:m])
    x .- x̅
end

# TODO: check if @fastmath really has an advantage
"""
    time_align(wf::AbstractVector{T}; 
        p::T=0.5, window::Tuple{Int, Int}=(500, 500), l::Int=300
    ) where {T <: AbstractFloat}

Find a window of size `sum(window) + 1` for the waveform `wf` such that 
the point where the `p`-th fraction of the maximum value is reached is 
at the center of the window. The maximum value is determined by the last 
`l` values.
"""
@fastmath function time_align(wf::AbstractVector{T}; p::T=0.5, 
window::Tuple{Int, Int}=(500, 500), l::Int=300) where {T<:AbstractFloat}
    @inbounds begin
        wfmax = mean(wf[end-l:end])
        idx = Int(find_crossing(smooth20(wf), p*wfmax, interpolate = false))
        if idx ≤ window[1] || length(wf) - window[2] < idx
            return wf
        end
        wf[idx-window[1]:idx+window[2]]
    end
end

@fastmath function decay_correction(wf::AbstractVector{T}, tau::T
) where {T<:AbstractFloat}
    rp = Vector{T}(undef, length(wf))
    rp[1] = wf[1]
    @inbounds for i in 2:length(rp)
        rp[i] = wf[i] + rp[i-1] - wf[i-1] * tau
    end
    rp
end

# TODO: replace by appropriate filters from DSP/RadiationDetectorDSP
function smooth20(wv::Vector{T}) where {T <: AbstractFloat}
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

@fastmath function find_crossing(wf::AbstractVector{T}, thold::T; 
window::Tuple{Int, Int} = (100, 1100), interpolate::Bool = true
)::T where {T<:AbstractFloat}
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

function chi_sq_wfs(x::AbstractVector{T}, y::AbstractVector{T}; 
window=400:600, l=150) where {T}
    varx = var(x[1:l])
    vary = var(y[1:l])
    sum((x[window] - y[window]).^2) / ((length(window) - 1)*(varx + vary))
end