# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).


using Plots # TODO: Temporary, replace plot commands by recipes!!!

function polaris_stdplots(df::DataFrame, idxs::Union{AbstractVector{<:Integer},Colon} = :)
    plot(
        histogram2d(
            mean.(df[idxs,:hit_x]), mean.(df[idxs,:hit_y]),
            bins = (-22:0.25:22, -22:0.25:22),
            xlabel = "x / mm", ylabel = "y / mm",
            ratio = 1,
            color = :viridis
        ),
        histogram2d(
            mean.(df[idxs,:hit_z]), mean.(df[idxs,:hit_y]),
            bins = (0:0.05:10, -22:0.25:22),
            xlabel = "z / mm", ylabel = "y / mm",
            ratio = 1,
            color = :viridis
        ),
        histogram2d(
            mean.(df[idxs,:hit_x]), mean.(df[idxs,:hit_z]),
            bins = (-22:0.25:22, 0:0.05:10),
            xlabel = "x / mm", ylabel = "z / mm",
            ratio = 1,
            color = :viridis
        ),
        stephist(
            sum.(df[idxs, :hit_edep]),
            bins = 0:2:700,
            xlabel = "E / keV", ylabel = "Counts"
        )
    )
end

polaris_stdplots(predicate::Function, ds::DataFrame, colname::Symbol...) =
    polaris_stdplots(ds, multifind(predicate, ds, colname...))

export polaris_stdplots
