# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

const x_CZT, y_CZT, z_CZT = 1000u"μm" .* (-4.5, 81.764-1.55, -62.8 + 19.75)
const δ = [x_CZT, y_CZT, z_CZT]
const α2 = 45.3695 * π/180
const δ2 = 1000u"μm" .* [64.1275058422767, -31.859041946909453, 0.4]


function get_global_cam_positions(c::NamedTuple)
    #transform local CZT coordinates to global coordinate system
    x_global = (cam = c.hit_x,)
    y_global = (cam = c.hit_y,)
    z_global = (cam = c.hit_z,)
    return x_global, y_global, z_global
end


# lowest allocation method
function transform_czt2_coords!(X::Vector{T}, Y::Vector{T}, Z::Vector{T}, 
motor_z::U)::Nothing where {T, U}
    t = SVector{3}(1000u"μm" .* (U(0), U(0), motor_z) .+ δ .+ δ2)
    R = RotZ(-α2)*RotY(π)
    @inbounds for i=eachindex(X)
        v = SVector{3, eltype(t)}(X[i], -Z[i], -Y[i])
        X[i], Y[i], Z[i] = round.(T, R * v + t)
    end
    nothing
end

# lowest allocation method 
function transform_czt1_coords!(X::Vector{T}, Y::Vector{T}, Z::Vector{T}, 
motor_z::U)::Nothing where {T, U}
    t = round.(T, (1000u"μm" .* (U(0), U(0), motor_z) .+ δ))
    @inbounds for i=eachindex(X)
        Y_i = Y[i]
        X[i] = t[1] + X[i]
        Y[i] = t[2] - Z[i]
        Z[i] = t[3] - Y_i
    end
    nothing
end

