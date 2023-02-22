mutable struct Cone{T <: AbstractFloat, TT <: QuantityMM}
    origin::Vector{TT} # in mm
    axis::Vector{T}
    α::T
    r1::Vector{T}
    r2::Vector{T}

    function Cone(::TT) where {T <: AbstractFloat, TT <: QuantityMM{T}}
        new{T,TT}(zeros(TT,3),zeros(T,3),0.0,zeros(T,3),zeros(T,3))
    end
end

#=
function Cone(origin::AbstractVector, axis::AbstractVector, α::Number, r1::AbstractVector, r2::AbstractVector)
    c = Cone(origin[1])
    c.origin = origin
    c.axis = axis/norm(axis)
    c.α = α
    c.r1 = r1/norm(r1)
    c.r2 = r2/norm(r2)
    return c
end
=#

function Cone(origin::Vector{QuantityMM{T}}, axis::Vector{QuantityMM{T}}, α::T) where {T <: AbstractFloat}
    c = Cone(origin[1])
    c.origin = origin
    c.axis = axis/norm(axis)
    c.α = α
    r1 = cross(T[0,1,0],c.axis)
    c.r1 = r1/norm(r1)
    r2 = cross(c.axis,c.r1)
    c.r2 = r2/norm(r2)
    return c
end
