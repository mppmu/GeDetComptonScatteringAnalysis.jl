# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

const regR = "R_(?<R>\\d+?\\.\\d+)mm"
const regPhi = "_Phi_(?<Phi>\\d+?\\.\\d+)deg"
const regT = "_T_(?<T>\\d+?\\.\\d+)K"
const regZ = "_Z_(?<Z>\\d+?\\.\\d+)mm"
const regM = "_measuretime_(?<time>\\d+?)sec"
const regV = "HV_(?<voltage>[\\d\\.]+?)V"

@inline getR(s::AbstractString)   = parse(Float64, match(Regex(regR), s).captures[1])
@inline getPhi(s::AbstractString) = parse(Float64, match(Regex(regPhi), s).captures[1])
@inline getT(s::AbstractString)   = parse(Float64, match(Regex(regT), s).captures[1])
@inline getZ(s::AbstractString)   = parse(Float64, match(Regex(regZ), s).captures[1])
@inline getM(s::AbstractString)   = parse(Float64, match(Regex(regM), s).captures[1])
@inline getV(s::AbstractString)   = parse(Float64, match(Regex(regV), s).captures[1])

############################
# filtered -> preprocessed #
############################
function fetch_relevant_filtered_files(
        sourcedir::AbstractString, phi::Number, z::Number, r::Number, 
        hv_in_fname::Bool, n_max_files::Int
    )::Vector{String}

    allfiles::Vector{String} = filter(file -> endswith(file, "filtered.h5")
                            && occursin("Phi_$(phi)", file)
                            && occursin("Z_$(z)", file)
                            && occursin("R_$(r)", file), readdir(sourcedir))
    if hv_in_fname
        allfiles = filter(file -> occursin("_HV_$(hv)", file), allfiles)
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

function build_output_prepocessed_file_name(files::Vector{String}, destdir::AbstractString, n::Int)::String
    destdir = joinpath(destdir, "")
    fpos, ftime = split(basename(files[end]), "measuretime_")
    fmtime, fdate = split(ftime, "sec")
    ending, _ = split(fdate, "filtered")
    fmtime = n*parse(Int32, fmtime)
    joinpath(destdir, 
        fpos * "measuretime_$(fmtime)sec" * ending * "preprocessed.lh5")
end

function write_preprocessed_file(f::LHDataStore, name::AbstractString , data::Tuple{Table, Table})::Nothing
    LegendHDF5IO.writedata(f.data_store, "$name", data[1])
    LegendHDF5IO.writedata(f.data_store, "czt", data[2])
    nothing
end


############################
# preprocessed -> results  #
############################
function is_preprocessed_file(f::AbstractString)::Bool 
    return (endswith(f, "-preprocessed.lh5") || throw(ArgumentError, "Expected prepocessed file as input.")) && isfile(f)
end

function read_preprocessed_file(f::AbstractString, name::AbstractString)::Tuple{Table, Table}
    @assert is_preprocessed_file(f)
    LHDataStore(f) do lhd
        lhd[name][:], lhd["czt"][:]
    end
end