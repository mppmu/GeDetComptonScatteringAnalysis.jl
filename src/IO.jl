# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

const regR = r"(?<=R_)\d+(?:\.\d+){0,1}mm(?=_)"
const regPhi = r"(?<=Phi_)\d+(?:\.\d+){0,1}deg(?=_)"
const regT = r"(?<=T_)\d+(?:\.\d+){0,1}K(?=_)"
const regZ = r"(?<=Z_)\d+(?:\.\d+){0,1}mm(?=_)"
const regM = r"(?<=measuretime_)\d+(?:\.\d+){0,1}s(?=ec)"
const regV = r"(?<=HV_)\d+(?:\.\d+){0,1}V(?=)"

@inline _parse(m::RegexMatch) = Float64(uparse(replace(m.match, "deg" => "°")))
@inline _parse(::Nothing) = NaN
@inline _match(r::Regex, s::AbstractString) = _parse(match(r, s))

@inline getR(s::AbstractString)   = _match(regR, s)
@inline getPhi(s::AbstractString) = _match(regPhi, s)
@inline getT(s::AbstractString)   = _match(regT, s)
@inline getZ(s::AbstractString)   = _match(regZ, s)
@inline getM(s::AbstractString)   = _match(regM, s)
@inline getV(s::AbstractString)   = _match(regV, s)

@inline compute_measuretime(files::Vector{String}) =
    Int32(ustrip(sum(map(getM, files))))

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
    return (endswith(f, "-filtered.h5") || throw(ArgumentError("Expected filtered file as input."))) && 
        (isfile(f) || throw(ArgumentError("File does not exist (missing path to file?)")))
end

function read_filtered_file(f::AbstractString, name::AbstractString)
    @assert is_filtered_file(f)
    LHDataStore(f) do lhd
        lhd[name][:], lhd["czt"][:], haskey(lhd, "czt2") ? lhd["czt2"][:] : missing
    end
end

function build_preprocessed_file_name(files::Vector{String}, 
        destdir::AbstractString)::String
    length(files) == 0 && throw(ArgumentError("No files were submitted"))
    fmtime = compute_measuretime(files)
    destdir = joinpath(destdir, "")
    fpos, ftime = split(basename(files[end]), "measuretime_")
    _, fdate = split(ftime, "sec")
    ending, _ = split(fdate, "filtered")
    joinpath(destdir, 
        fpos * "measuretime_$(fmtime)sec" * ending * "preprocessed.lh5")
end

function write_preprocessed_file(
        fileout::AbstractString, det_name::AbstractString,
        data::Tuple{detTable, cztTable, Int, typeof(Cs_energy)}, 
        mode = "cw")::Nothing
    LHDataStore(fileout, mode) do f
        LegendHDF5IO.writedata(f.data_store, "$det_name", data[1])
        LegendHDF5IO.writedata(f.data_store, "czt", data[2])
        f["idx_c"] = data[3]
        f["econv"] = data[4]
    end
    nothing
end


############################
# preprocessed -> results  #
############################
function is_preprocessed_file(f::AbstractString)::Bool 
    return (endswith(f, "-preprocessed.lh5") || throw(ArgumentError("Expected prepocessed file as input."))) && 
        (isfile(f) || throw(ArgumentError("File does not exist (missing path to file?)")))
end

function read_preprocessed_file(f::AbstractString, name::AbstractString)::Tuple{detTable, cztTable, Int, typeof(Cs_energy)}
    @assert is_preprocessed_file(f)
    LHDataStore(f) do lhd
        lhd[name][:], lhd["czt"][:], lhd["idx_c"], lhd["econv"]
    end
end