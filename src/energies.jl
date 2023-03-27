const Cs_energy = 661.66u"keV"

# const econv = Dict{Float64, Float64}( 
#     300.0 => 0.000746364303766692,
#     # 600.0 => 0.0007730410437211994,
#     600.0 => 0.0006952481325306902,
#     # 900.0 => 0.0006666592487136949,
#     900.0 => 0.0006781005477809491,
#     1200.0 => 0.0006189596091123743, 
#     3500.0 => 0.0005490928284752016,
#     3000.0 => 0.000653
# )


# ----- INFO ------
# Energies should be determined from pulses...
function correct_DAQ_energy!(det::detTable, cf::T)::Nothing where {T <: AbstractFloat}
    for i in eachindex(det)
        new_e::T = DAQ_energy_corr(det.samples[i], cf * det.DAQ_energy[i])
        det.DAQ_energy[i] = new_e / cf
    end
    nothing
end

function DAQ_energy_corr(wf::Vector{UInt16}, daqe::T)::T where {T <: AbstractFloat}
    fit_par::Vector{T} = T[0.3762989397047389, 0.7729112545305099, 0.0007296045244350602]
    bs::T = baseline_slope(wf)
    if bs < T(-1.5)
        return daqe - par(bs, fit_par)
    else
        return daqe
    end
end   

baseline_slope(wf::Vector{UInt16}) = mean(view(wf, 1400:1500)) - mean(view(wf, 1:100))
par(x::T, p::Vector{T}) where {T <: AbstractFloat} = p[1] + p[2] * x + p[3] * x^2


"""
    get_econv(sourcedir [; idx_c=1 nseg=5, max=1_000_000, name="segBEGe"])

compute energy conversion factor. From all files using only the 
DAQ_energy from the core. `idx_c` denotes the index of the core segment, 
`nseg` the number of segments ,`max` the maximal raw value to 
consider in the histogram and `name` the name of the table
"""
function get_econv(det::detTable; idx_c::Int = 1, bsize::Int = 1000, 
        max::Int = 1_500_000, verbose::Bool = true)::typeof(Cs_energy)
    h::Histogram{Float64} = fit(
        Histogram{Float64}, 
        det.DAQ_energy[findall(det.chid .== idx_c)], 
        (1:bsize:max))
    _, peakpos = RadiationSpectra.peakfinder(h)
    cf = Cs_energy / maximum(peakpos)
    verbose && @info "energy conversion factor" cf
    return cf
end

"""
    compute_measured_pulse_amplitudes(det::detTable, 
    baseline_length::RealQuantity, tail_start::RealQuantity, 
    tail_stop::RealQuantity, tau_window::AbstractRange{U}, 
    trap_params::NTuple{2, T}) where {T<:RealQuantity, U<:RealQuantity}


"""
function compute_measured_pulse_amplitudes(det, 
baseline_length::RealQuantity, tail_start::RealQuantity, 
tail_stop::RealQuantity, tau_window::AbstractRange{U}, 
trap_params::NTuple{2, T}) where {T<:RealQuantity, U<:RealQuantity}
    # samples = convert.(Vector{Float64}, det.samples)
    samples = det.samples
    taus_and_means = StructArray(
        get_tau.(samples, baseline_length, tail_start, tail_stop))
    T1 = eltype(taus_and_means.tau)
    # determine tau from peak/gaussian/cauchy fit
    h = fit(Histogram, ustrip(taus_and_means.tau), tau_window)
    # simple peak fit
    _, tau_peaks = RadiationSpectra.peakfinder(h)
    tau = T1(tau_peaks[end])
    # Cauchy fit 
    # param = curve_fit(cauchy_model, tau_window[2:end], h.weights, p0)
    # tau = param[2]

    # define relevant filters
    deconvflt = InvCRFilter(tau)
    trapflt = TrapezoidalChargeFilter(trap_params...)
    samples = trapflt.(deconvflt.(shift_waveform.(
        samples, -taus_and_means.baseline_mean)))
    findmax.(samples)
end