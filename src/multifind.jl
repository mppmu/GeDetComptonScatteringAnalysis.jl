# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).


function multifind(predicate::Function, ds::DataFrame, colname::Symbol...)
    columns = map(c -> ds[c], colname)
    f(xs) = predicate(xs...)
    find(f, zip(columns...))
end
