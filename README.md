# GeDetComptonScatteringAnalysis.jl

## Documentation

Main functions of interest are `stack_and_merge_at_z()` and `get_all_z()`.
`stack_and_merge_at_z` produces for a specified radius, height z and 
angle phi one filtered lh5 file, containing the core DAQ_energies and 
merged data from both czt cameras.\
`get_all_z` uses those previously produced files to reconstruct the z 
position of each 'two-hit' event.

## Examples

Produce filtered files

```julia
julia> datapath = "path/to/data/files/"
julia> destdir = "path/to/directory/of/filtered/files/"
julia> hv = 600.        # applied Voltage
julia> r = 47.8         # radial position of radiation source (in mm)
julia> phi = 88.5       # azimuthal position of radiation source (in degrees)
julia> z = 62.8         # height of czt cameras (in mm)
julia> name = "segBEGe" # group name of detector in raw lh5 files
julia> stack_and_merge_at_z(datapath, destdir, r, phi, z, hv, name)
```
reconstruct z positions
```julia
julia> mtime, R, z = get_all_z(destdir)
```
`mtime` contains the effective measuretime for all radii `R` with corresponding
reconstructed `z` from two hit events.