# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

const x_CZT, y_CZT, z_CZT = 1000u"μm" .* (-4.5, 81.764-1.55, -62.8 + 19.75)
const δ = [x_CZT, y_CZT, z_CZT]
const α2 = 45.3695 * π/180
const δ2 = 1000u"μm" .* [64.1275058422767, -31.859041946909453, 0.4]
const cntr = 81.76361317572471u"mm"

# lowest allocation method
function transform_czt2_coords!(X::Vector{T}, Y::Vector{T}, Z::Vector{T}, 
motor_z::TT)::Nothing where {T, TT <: QuantityMM{Float64}}
    t = SVector{3}((zero(T), zero(T), round(T,motor_z)) .+ δ .+ δ2)
    R = RotZ(-α2)*RotY(π)
    @inbounds for i=eachindex(X)
        v = SVector{3, eltype(t)}(X[i], -Z[i], -Y[i])
        X[i], Y[i], Z[i] = round.(T, R * v + t)
    end
    nothing
end

# lowest allocation method 
function transform_czt1_coords!(X::Vector{T}, Y::Vector{T}, Z::Vector{T}, 
motor_z::TT)::Nothing where {T, TT <: QuantityMM{Float64}}
    t = round.(T, (zero(T), zero(T), motor_z) .+ δ)
    @inbounds for i=eachindex(X)
        Y_i = Y[i]
        X[i] = t[1] + X[i]
        Y[i] = t[2] - Z[i]
        Z[i] = t[3] - Y_i
    end
    nothing
end

function merge_cameras_and_transform_coordinates(::AbstractDetTable, czt::cztTable, ::Missing, motor_z::QuantityMM{Float64})::cztTable
    transform_czt1_coords!(czt.hit_x.data, czt.hit_y.data, czt.hit_z.data, motor_z)
    czt
end

function merge_cameras_and_transform_coordinates(det::AbstractDetTable, czt::cztTable, czt2::cztTable, motor_z::QuantityMM{Float64})::cztTable
    # TODO: maybe consider doing transformation directly in main loop here
    # instead of seperately for each camera
    transform_czt1_coords!(czt.hit_x.data, czt.hit_y.data, czt.hit_z.data, motor_z)
    transform_czt2_coords!(czt2.hit_x.data, czt2.hit_y.data, czt2.hit_z.data, motor_z)
    N::Int, N1::Int, N2::Int = length(det), length(czt), length(czt2)
    total_cam_hits::Int = length(czt.hit_x.data) + length(czt2.hit_x.data)

    # simple Vectors
    evt_nhits = Vector{Int32}(undef, N)
    evt_t = Vector{eltype(czt.evt_t)}(undef, N)
    evt_issync::BitVector = falses(N)

    # data for VectorOfVectors
    elem_ptr = Vector{Int64}(undef, N+1)
    detno = Vector{Int32}(undef, total_cam_hits)
    edep = Vector{eltype(czt.hit_edep.data)}(undef, total_cam_hits)
    t = Vector{eltype(czt.hit_t.data)}(undef, total_cam_hits)
    x = Vector{eltype(czt.hit_x.data)}(undef, total_cam_hits)
    y = Vector{eltype(czt.hit_y.data)}(undef, total_cam_hits)
    z = Vector{eltype(czt.hit_z.data)}(undef, total_cam_hits)

    i1, i2, ep = 1, 1, 0
    for i=eachindex(det)
        elem_ptr[i] = ep + 1
        evt_no_cam1, evt_no_cam2 = czt.evt_no[i1], czt2.evt_no[i2]
        # events with hits in both cameras
        if evt_no_cam1 == evt_no_cam2 == det.evt_no[i]
            evt_nhits[i] = czt.evt_nhits[i1] + czt2.evt_nhits[i2]
            evt_t[i] = min(czt.evt_t[i1], czt2.evt_t[i2])
            for n=Base.OneTo(czt.evt_nhits[i1])
                detno[ep + n] = czt.hit_detno[i1][n]
                edep[ep + n] = czt.hit_edep[i1][n]
                t[ep + n] = czt.hit_t[i1][n]
                x[ep + n] = czt.hit_x[i1][n]
                y[ep + n] = czt.hit_y[i1][n]
                z[ep + n] = czt.hit_z[i1][n]
            end
            ep += czt.evt_nhits[i1]
            for n=Base.OneTo(czt2.evt_nhits[i2])
                detno[ep + n] = 4 + czt2.hit_detno[i2][n]
                edep[ep + n] = czt2.hit_edep[i2][n]
                t[ep + n] = czt2.hit_t[i2][n]
                x[ep + n] = czt2.hit_x[i2][n]
                y[ep + n] = czt2.hit_y[i2][n]
                z[ep + n] = czt2.hit_z[i2][n]
            end
            ep += czt2.evt_nhits[i2]
            i1 += 1
            i2 += 1
        # events with hits in only camera 1
        elseif evt_no_cam1 == det.evt_no[i]
            evt_nhits[i] = czt.evt_nhits[i1]
            evt_t[i] = czt.evt_t[i1]
            for n=Base.OneTo(czt.evt_nhits[i1])
                detno[ep + n] = czt.hit_detno[i1][n]
                edep[ep + n] = czt.hit_edep[i1][n]
                t[ep + n] = czt.hit_t[i1][n]
                x[ep + n] = czt.hit_x[i1][n]
                y[ep + n] = czt.hit_y[i1][n]
                z[ep + n] = czt.hit_z[i1][n]
            end
            ep += czt.evt_nhits[i1]
            i1 += 1
        # events with hits in only camera 2
        elseif evt_no_cam2 == det.evt_no[i]
            evt_nhits[i] = czt2.evt_nhits[i2]
            evt_t[i] = czt2.evt_t[i2]
            for n=Base.OneTo(czt2.evt_nhits[i2])
                detno[ep + n] = 4 + czt2.hit_detno[i2][n]
                edep[ep + n] = czt2.hit_edep[i2][n]
                t[ep + n] = czt2.hit_t[i2][n]
                x[ep + n] = czt2.hit_x[i2][n]
                y[ep + n] = czt2.hit_y[i2][n]
                z[ep + n] = czt2.hit_z[i2][n]
            end
            ep += czt2.evt_nhits[i2]
            i2 += 1
        # events with no hits in the cameras
        # this should not happen, as the filtered files should
        # only contain events with hits in the detector AND the camera
        else 
            @warn "Event with event number $(det.evt_no[i]) in the detector has no counterpart in the cameras."
        end
        i1 > N1 && (i1 = N1)            # TODO: maybe get rid of these?
        i2 > N2 && (i2 = N2)
    end
    elem_ptr[end] = total_cam_hits + 1
    Table(
        evt_no = det.evt_no,
        evt_t = evt_t,
        evt_nhits = evt_nhits,
        evt_issync = evt_issync,
        hit_edep = VectorOfVectors(edep, elem_ptr, no_consistency_checks),
        hit_t = VectorOfVectors(t, elem_ptr, no_consistency_checks),
        hit_detno = VectorOfVectors(detno, elem_ptr, no_consistency_checks),
        hit_x = VectorOfVectors(x, elem_ptr, no_consistency_checks),
        hit_y = VectorOfVectors(y, elem_ptr, no_consistency_checks),
        hit_z = VectorOfVectors(z, elem_ptr, no_consistency_checks),
    )[:]
end
