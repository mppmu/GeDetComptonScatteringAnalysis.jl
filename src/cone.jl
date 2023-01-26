measuredd = 106.5
measuredh = 1.3
z_offset = -25.34
global campos = [0.0,measuredd+z_offset,66.4+measuredh-47.5]u"mm"
export campos

mutable struct Cone
    origin::AbstractVector
    axis::AbstractVector
    α::Number
    r1::AbstractVector
    r2::AbstractVector

    function Cone()
        new([0,0,0],[0,0,0],0,[0,0,0],[0,0,0])
    end
end

function Cone(origin::AbstractVector, axis::AbstractVector, α::Number, r1::AbstractVector, r2::AbstractVector)
    c = Cone()
    c.origin = origin
    c.axis = axis/norm(axis)
    c.α = α
    c.r1 = r1/norm(r1)
    c.r2 = r2/norm(r2)
    return c
end

function Cone(origin::AbstractVector, axis::AbstractVector, α::Number)
    c = Cone()
    c.origin = origin
    c.axis = axis/norm(axis)
    c.α = α
    r1 = cross(Vector([0,1,0]),c.axis)
    c.r1 = r1/norm(r1)
    r2 = cross(c.axis,c.r1)
    c.r2 = r2/norm(r2)
    return c
end

export Cone


#=
function cone_points(cone::Cone,t::Number,h::Number)::AbstractVector #origin,t,h,H,r1,r2,α)::Vector
    return cone.origin + cos(cone.α)*h*cone.axis + h*sin(cone.α)*(cone.r1*cos(t) + cone.r2*sin(t))
end

function plot_cone(cone::Cone, p::AbstractVector)

    d = norm(cone.origin - p)

    Cx = [cone_points(cone,t,h)[1] for t in 0:(2*pi/100):(2*pi+0.01) for h in [0,d]]
    Cy = [cone_points(cone,t,h)[2] for t in 0:(2*pi/100):(2*pi+0.01) for h in [0,d]]
    Cz = [cone_points(cone,t,h)[3] for t in 0:(2*pi/100):(2*pi+0.01) for h in [0,d]]

    Dx = [cone_points(cone,t,d)[1] for t in 0:(2*pi/100):(2*pi+0.01)]
    Dy = [cone_points(cone,t,d)[2] for t in 0:(2*pi/100):(2*pi+0.01)]
    Dz = [cone_points(cone,t,d)[3] for t in 0:(2*pi/100):(2*pi+0.01)]

    plot!(Cx, Cy, Cz, label = "", color = :gray, width = 0.1)
    plot!(Dx, Dy, Dz, label = "", color = :gray, width = 1.5)
    scatter!([p[1]],[p[2]],[p[3]], markersize = 10, color = :gray)

end

export plot_cone
=#
