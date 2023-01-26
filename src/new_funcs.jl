export stack_and_merge_at_z, get_all_z


"""
    get_all_z(sourcedir; center=CNTR_VALUE)

loop through all files in `sourcedir` and return for each file the radius, 
measuretime and reconstructed z's
"""
function get_all_z(sourcedir; name="segBEGe", center=81.76361317572471)
    ffiles = readdir(sourcedir)
    R = zeros(length(ffiles))
    mtime = zeros(length(ffiles))
    z = Vector{Float64}[]
    for i=eachindex(ffiles)
        R[i] = getR(ffiles[i])
        mtime[i] = getM(ffiles[i])
        rec_zs = getz(joinpath(sourcedir, ffiles[i]); name, center, ew=20)
        push!(z, rec_zs[:, 1])
    end
    mtime, R, z
end

