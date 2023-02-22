# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

using Plots # TODO: Temporary, replace plot commands by recipes!!!

function multifind(predicate::Function, data, colname::Symbol...)
    cols = Tables.columns(data)
    f(xs) = predicate(xs...)
    findall(f, collect(zip(cols...))) # TODO: Get rid of collect here
end

function polaris_stdplots(events, idxs::Union{AbstractVector{<:Integer},Colon} = :)
    hit_x = deepmap(x -> ustrip(uconvert(u"mm", x)), events[idxs].hit_x)
    hit_y = deepmap(x -> ustrip(uconvert(u"mm", x)), events[idxs].hit_y)
    hit_z = deepmap(x -> ustrip(uconvert(u"mm", x)), events[idxs].hit_z)
    hit_edep = deepmap(x -> ustrip(uconvert(u"keV", x)), events[idxs].hit_edep)

    plot(
        histogram2d(
            mean.(hit_x), mean.(hit_y),
            bins = (-22:0.25:22, -22:0.25:22),
            xlabel = "x / mm", ylabel = "y / mm",
            ratio = 1,
            color = :viridis
        ),
        histogram2d(
            mean.(hit_z), mean.(hit_y),
            bins = (0:0.05:10, -22:0.25:22),
            xlabel = "z / mm", ylabel = "y / mm",
            ratio = 1,
            color = :viridis
        ),
        histogram2d(
            mean.(hit_x), mean.(hit_z),
            bins = (-22:0.25:22, 0:0.05:10),
            xlabel = "x / mm", ylabel = "z / mm",
            ratio = 1,
            color = :viridis
        ),
        stephist(
            sum.(hit_edep),
            bins = 0:2:700,
            xlabel = "E / keV", ylabel = "Counts"
        )
    )
end

polaris_stdplots(predicate::Function, data, colname::Symbol...) =
    polaris_stdplots(data, multifind(predicate, data, colname...))

export polaris_stdplots

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
