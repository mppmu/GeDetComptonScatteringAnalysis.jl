# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).


function multifind(predicate::Function, data, colname::Symbol...)
    cols = Tables.columns(data)
    f(xs) = predicate(xs...)
    findall(f, collect(zip(cols...))) # TODO: Get rid of collect here
end
