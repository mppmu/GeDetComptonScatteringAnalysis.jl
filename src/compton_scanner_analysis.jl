in_mm(x::Number) = ustrip(uconvert(u"mm", x))

function get_z_from_2_hit_events(s::NamedTuple, c::NamedTuple, R::Number; Δz = 2)

    #get z coordinates for scatter point in segBEGe
    zθ = get_z_from_energies(s,c,R)
    zα = get_z_from_camera(s,c,R)

    #compare them: if they agree within Δz, keep them
    for z in zα
        if abs(z - zθ) < Δz
            return (true, zθ, z)
        end
    end

    #else discard them
    return (false, zθ)
end

export get_z_from_2_hit_events


function get_z_from_energies(s::NamedTuple, c::NamedTuple, R::Number; campos::AbstractVector = campos)
    try
        x_global, y_global, z_global = get_global_cam_positions(c, campos)
        θ = compton_angle(s.energy+sum(c.hit_edep), sum(c.hit_edep))
        return in_mm(z_global.cam[1] + sqrt(x_global.cam[1] ^ 2 + (y_global.cam[1] - R*u"mm") ^ 2) * cot(θ))
    catch e
        error(e)
    end
end

export get_z_from_energies


function get_z_from_camera(s::NamedTuple, c::NamedTuple, R::Number; campos::AbstractVector = campos)

    try
        x_global, y_global, z_global = get_global_cam_positions(c, campos)
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

    catch e
        error(e)
    end
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


function get_global_cam_positions(c::NamedTuple, campos::AbstractVector)
    #transform local CZT coordinates to global coordinate system
    x_global = (cam = (c.hit_x) .+ campos[1],)
    y_global = (cam = -1*(c.hit_z) .+ campos[2],)
    z_global = (cam = -1*(c.hit_y) .+ campos[3],)
    return x_global, y_global, z_global
end

function swap_CZT_hits(c::NamedTuple)
    return (hit_x = [c.hit_x[2], c.hit_x[1]], hit_y = [c.hit_y[2], c.hit_y[1]], hit_z = [c.hit_z[2], c.hit_z[1]], hit_edep = [c.hit_edep[2], c.hit_edep[1]])
end

export swap_CZT_hits
