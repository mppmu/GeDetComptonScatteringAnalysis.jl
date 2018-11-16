# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).


function multifind(predicate::Function, ds::DataFrame, colname::Symbol...)
    cols = map(c -> ds[c], colname)
    f(xs) = predicate(xs...)
    findall(f, collect(zip(cols...))) # TODO: Get rid of collect here
end
