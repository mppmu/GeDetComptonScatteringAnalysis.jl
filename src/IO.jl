# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

regR = "R_(?<R>\\d+?\\.\\d+)mm"
regPhi = "_Phi_(?<Phi>\\d+?\\.\\d+)deg"
regT = "_T_(?<T>\\d+?\\.\\d+)K"
regZ = "_Z_(?<Z>\\d+?\\.\\d+)mm"
regM = "_measuretime_(?<time>\\d+?)sec"
regV = "HV_(?<voltage>[\\d\\.]+?)V"
getR(s::String) = parse(Float64, match(Regex(regR), s).captures[1])
getPhi(s::String) = parse(Float64, match(Regex(regPhi), s).captures[1])
getT(s::String) = parse(Float64, match(Regex(regT), s).captures[1])
getZ(s::String) = parse(Float64, match(Regex(regZ), s).captures[1])
getM(s::String) = parse(Float64, match(Regex(regM), s).captures[1])
getV(s::String) = parse(Float64, match(Regex(regV), s).captures[1])

############################
# filtered -> preprocessed #
############################
function fetch_relevant_filtered_files(sourcedir, phi, z, r, hv_in_fname, n_max_files)
    allfiles = filter(file -> endswith(file, "filtered.h5")
                            && occursin("Phi_$(phi)", file)
                            && occursin("Z_$(z)", file)
                            && occursin("R_$(r)", file), readdir(sourcedir))
    if hv_in_fname
        allfiles = filter(file -> contains(file, "_HV_$hv"), allfiles)
    end
    if n_max_files > 0
        allfiles = allfiles[1:min(length(allfiles), n_max_files)]
    end
    joinpath.(sourcedir, allfiles)
end

function is_filtered_file(f::AbstractString)::Bool
    return (endswith(f, "-filtered.h5") || throw(ArgumentError, "Expected filtered file as input.")) && isfile(f)
end

function read_filtered_file(f::AbstractString, name::AbstractString)
    @assert is_filtered_file(f)
    LHDataStore(f) do lhd
        lhd[name][:], lhd["czt"][:], haskey(lhd, "czt2") ? lhd["czt2"][:] : missing
    end
end

function build_output_prepocessed_file_name(files, destdir, n)
    destdir = joinpath(destdir, "")
    fpos, ftime = split(basename(files[end]), "measuretime_")
    fmtime, fdate = split(ftime, "sec")
    ending, _ = split(fdate, "filtered")
    fmtime = n*parse(Int32, fmtime)
    joinpath(destdir, 
        fpos * "measuretime_$(fmtime)sec" * ending * "preprocessed.lh5")
end

write_preprocessed_file(f::LHDataStore, name::String , data::Tuple{Table, Table}) = begin
    LegendHDF5IO.writedata(f.data_store, "$name", data[1])
    LegendHDF5IO.writedata(f.data_store, "czt", data[2])
end


############################
# preprocessed -> results  #
############################
function is_preprocessed_file(f::AbstractString)::Bool 
    return (endswith(f, "-preprocessed.lh5") || throw(ArgumentError, "Expected prepocessed file as input.")) && isfile(f)
end

function read_preprocessed_file(f::AbstractString, name::AbstractString)
    @assert is_preprocessed_file(f)
    LHDataStore(f) do lhd
        lhd[name][:], lhd["czt"][:]
    end
end