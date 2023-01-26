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

