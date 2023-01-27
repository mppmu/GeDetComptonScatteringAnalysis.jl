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

@inline function is_valid_2hit(evt)
    length(evt.hit_x) == 2 && hypot(evt.hit_x[2] - evt.hit_x[1],
        evt.hit_y[2] - evt.hit_y[1],
        evt.hit_z[2] - evt.hit_z[1]) > 3.0u"mm"
end