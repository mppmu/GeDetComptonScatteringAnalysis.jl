# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

const regR = r"(?<=R_)\d+(?:\.\d+){0,1}(?=mm)"
const regPhi = r"(?<=Phi_)\d+(?:\.\d+){0,1}(?=deg)"
const regT = r"(?<=T_)\d+(?:\.\d+){0,1}(?=K)"
const regZ = r"(?<=Z_)\d+(?:\.\d+){0,1}(?=mm)"
const regM = r"(?<=measuretime_)\d+(?:\.\d+){0,1}(?=sec)"
const regV = r"(?<=HV_)\d+(?:\.\d+){0,1}(?=V)"

_parse(m::RegexMatch) = parse(Float64, m.match)
_parse(::Nothing) = NaN
_match(r::Regex, s::AbstractString)::Float64 = _parse(match(r, s))

@inline getR(s::AbstractString)   = _match(regR, s)
@inline getPhi(s::AbstractString) = _match(regPhi, s)
@inline getT(s::AbstractString)   = _match(regT, s)
@inline getZ(s::AbstractString)   = _match(regZ, s)
@inline getM(s::AbstractString)   = _match(regM, s)
@inline getV(s::AbstractString)   = _match(regV, s)

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
    return (endswith(f, "-filtered.h5") || throw(ArgumentError, "Expected filtered file as input.")) && 
           (isfile(f) || throw(ArgumentError, "File does not exist (missing path to file?)"))
end

function read_filtered_file(f::AbstractString, name::AbstractString)
    @assert is_filtered_file(f)
    LHDataStore(f) do lhd
        lhd[name][:], lhd["czt"][:], haskey(lhd, "czt2") ? lhd["czt2"][:] : missing
    end
end

function build_preprocessed_file_name(files::Vector{String}, destdir::AbstractString, n::Int)::String
    destdir = joinpath(destdir, "")
    fpos, ftime = split(basename(files[end]), "measuretime_")
    fmtime, fdate = split(ftime, "sec")
    ending, _ = split(fdate, "filtered")
    fmtime = n*parse(Int32, fmtime)
    joinpath(destdir, 
        fpos * "measuretime_$(fmtime)sec" * ending * "preprocessed.lh5")
end

function write_preprocessed_file(f::LHDataStore, name::AbstractString , data::Tuple{detTable, cztTable})::Nothing
    LegendHDF5IO.writedata(f.data_store, "$name", data[1])
    LegendHDF5IO.writedata(f.data_store, "czt", data[2])
    nothing
end


############################
# preprocessed -> results  #
############################
function is_preprocessed_file(f::AbstractString)::Bool 
    return (endswith(f, "-preprocessed.lh5") || throw(ArgumentError, "Expected prepocessed file as input.")) &&
           (isfile(f) || throw(ArgumentError, "File does not exist (missing path to file?)"))
end

function read_preprocessed_file(f::AbstractString, name::AbstractString)::Tuple{detTable, cztTable}
    @assert is_preprocessed_file(f)
    LHDataStore(f) do lhd
        lhd[name][:], lhd["czt"][:]
    end
end