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