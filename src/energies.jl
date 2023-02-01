const Cs_energy = 661.66u"keV"

const econv = Dict( 
                    300.0 => 0.000746364303766692,
                    # 600.0 => 0.0007730410437211994,
                    600.0 => 0.0006952481325306902,
                    # 900.0 => 0.0006666592487136949,
                    900.0 => 0.0006781005477809491,
                    1200.0 => 0.0006189596091123743, 
                    3500.0 => 0.0005490928284752016,
                    3000.0 => 0.000653
                )


# ----- INFO ------
# Energies should be determined from pulses...
function correct_DAQ_energy!(det, cf)
    for i=eachindex(det)
        new_e = DAQ_energy_corr(det.samples[i], cf*det.DAQ_energy[i])
        det.DAQ_energy[i] = new_e / cf
    end
end

function DAQ_energy_corr(wf, daqe)
    fit_par = [0.3762989397047389, 0.7729112545305099, 0.0007296045244350602]
    bs = baseline_slope(wf)
    if bs < -1.5
        return daqe - par(bs, fit_par)
    else
        return daqe
    end
end   

@. par(x, p) = p[1] + p[2]x + p[3]x^2


"""
    get_econv(sourcedir [;i0=1 nseg=5, max=1_000_000, name="segBEGe"])

compute energy conversion factor. From all files using only the 
DAQ_energy from the core. `i0` denotes the index of the core segment, 
`nseg` the number of segments ,`max` the maximal raw value to 
consider in the histogram and `name` the name of the table
"""
function get_econv(sourcedir; i0=1, bsize=1000, max=1_500_000, name="segBEGe")::Float64
    E = Int32[]
    files = readdir(sourcedir)
    N = length(files)
    for i=Base.OneTo(N)
        print("\e[0;0H\e[2J")
        println("file $i/$N")
        try
            lhd = LHDataStore(joinpath(sourcedir, files[i]))
            core_e = get_daqe(lhd[name], i0)
            append!(E, core_e)
        catch
            nothing
        end
    end
    h = fit(Histogram{Float64}, E, (1:bsize:max))
    _, peakpos = RadiationSpectra.peakfinder(h)
    return 661.660 / peakpos[1]
end

function get_daqe(x::TypedTables.Table, i::Int)
    daqe = x.DAQ_energy[:]
    chid = x.chid[:]
    daqe[findall(x -> x == i, chid)]
end

export get_daqe