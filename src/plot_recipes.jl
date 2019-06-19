# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).


using Plots # TODO: Temporary, replace plot commands by recipes!!!

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
