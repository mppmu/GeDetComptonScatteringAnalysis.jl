# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

function get_z_from_2_hit_events(
        det::detTuple, czt::cztTuple, R::QuantityMM{Float64}; 
        Δz::QuantityMM = 2.0u"mm", hv::Float64
    )::Tuple{Int, QuantityMM{Float64}, QuantityMM{Float64}}

    #if econv*det.DAQ_energy < 250
    #    return (0, NaN*u"mm", NaN*u"mm")
    #end

    #get z coordinates for scatter point in detector
    #compare them: if they agree within Δz, keep them
    zθ::QuantityMM{Float64} = get_z_from_energies(det, czt, R, hv)
    zα::QuantityMM{Float64} = get_z_from_camera(czt, R)
    if !isnan(zα) && abs(zα - zθ) < Δz 
        return (1, zθ, zα)
    end
    
    #swap hits and try again
    czts = swap_CZT_hits(czt)
    zθ = get_z_from_energies(det, czts, R, hv)
    zα = get_z_from_camera(czts, R)
    if !isnan(zα) && abs(zα - zθ) < Δz 
        return (2, zθ, zα)
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


function get_z_from_camera(czt::cztTuple, R::QuantityMM{Float64})::QuantityMM{Float64}

    let x = czt.hit_x, y = czt.hit_y, z = czt.hit_z
        α::Float64 = compton_angle(sum(czt.hit_edep), czt.hit_edep[2])
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
                    return z1
                elseif validate_z(z2,cone,R)
                    return z2
                end
            catch e
                #this catches all imaginary solutions in get_possible_z_from_camera
                if !(e isa DomainError) error(e) end
            end
        end
        return NaN * u"mm"
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