# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).
using RadiationDetectorDSP: SamplesOrWaveform

function triangular_dither(x::Number, width::Number = one(typeof(x)) * unit(x))
    T = float(typeof(ustrip(x)))
    r = rand(T)
    tr = (r >= T(0.5)) ? - sqrt(2 - 2*r) + 1 : sqrt(2*r) - 1
    x + tr * width
end


function polaris_dither(data)
    T_time = typeof(Float64(1)u"μs")
    T_pos = typeof(Float32(1)u"mm")
    T_energy = typeof(Float64(1)u"keV")

    cols = Tables.columns(data)

    Tables.materializer(data)(merge(cols, (
        evt_t = convert.(T_time, data.evt_t),
        hit_x = deepmap(x -> convert(T_pos, triangular_dither(x)), data.hit_x),
        hit_y = deepmap(x -> convert(T_pos, triangular_dither(x)), data.hit_y),
        hit_z = deepmap(x -> convert(T_pos, triangular_dither(x)), data.hit_z),
        hit_edep = deepmap(x -> convert(T_energy, triangular_dither(x)), data.hit_edep),
    )))
end


function correct_timestamps!(tables...)
    time_unit = unit(eltype(tables[1].evt_t))
    t_sync = map(tbl -> ustrip.(uconvert.(time_unit, tbl.evt_t[findall(tbl.evt_issync)])), tables)
    n_sync_events = min(map(length, t_sync)...)
    resize!.(t_sync, n_sync_events)

    all(x -> length(t_sync) - n_sync_events <= 2, t_sync) || @error "Number of sync events doesn't match"

    t_offs = tables[1].evt_t[1]

    for i in 2:length(tables)
        # Determine timestamp mapping relative to reference table:
        t_sync_fit = curve_fit((x, p) -> p[1] .+ p[2] * x, t_sync[i], t_sync[1], [1.0, 1.0])
        sync_p1, sync_p2 = t_sync_fit.param[1], t_sync_fit.param[2]

        # Correct timestamps:
        tables[i].evt_t .= sync_p1 * time_unit .+ sync_p2 .* tables[i].evt_t .- t_offs
    end

    tables[1].evt_t .= tables[1].evt_t .- t_offs

    tables
end


function find_common_events(tables::Tuple{Any,Any}, delta_t::Number)
    t = (tables[1].evt_t, tables[2].evt_t)
    sel = Vector{Int}(), Vector{Int}()

    idxs = (eachindex(t[1]), eachindex(t[2]))
    i = map(firstindex, idxs)
    i_max = map(lastindex, idxs)

    while i[1] < i_max[1] && i[2] < i_max[2]
        t1 = t[1][i[1]]
        t2 = t[2][i[2]]
        if abs(t1 - t2) < delta_t
            push!(sel[1], i[1])
            push!(sel[2], i[2])
            i = (i[1] + 1, i[2] + 1)
        elseif t1 <= t2
            i = (i[1] + 1, i[2])
        else
            i = (i[1], i[2] + 1)
        end
    end
    sel
end

const electron_mass = 510.9989u"keV"

function compton_angle(E_in::Number, E_out::Number)::Float64
    cos_theta::Float64 = 1 - uconvert(Unitful.NoUnits, electron_mass * (1/E_out - 1/E_in))
    (-1 < cos_theta < 1) ? acos(cos_theta) : NaN
end

compton_E_out(E_in::Number, θ::Real) = E_in / (1 + E_in/electron_mass * (1 - cos(θ)))

@fastmath function linear_regression(x::Vector{<:Real}, y::Vector{<:Real})::Tuple # Substitutes linear fit --> much faster
    @assert length(x) == length(y) "x and y must have the same length."
    T=Float64
    x_mean::T = mean(x)
    y_mean::T = mean(y)
    num::T = 0.0
    nom::T = 0.0
    for i in eachindex(x)
        x_res = (x[i] - x_mean)
        num += x_res * (y[i] - y_mean)
        nom += x_res * x_res
    end
    slope::T = num / nom
    offset::T = y_mean - slope * x_mean
    return offset, slope
end

@fastmath function linear_regression(x::AbstractVector{T}, y::AbstractVector{T})::Tuple{T, T} where {T <: AbstractFloat} # Substitutes linear fit --> much faster
    @assert length(x) == length(y) "x and y must have the same length."
    x_mean::T = mean(x)
    y_mean::T = mean(y)
    num::T = 0.0
    nom::T = 0.0
    @inbounds for i in eachindex(x)
        x_res = (x[i] - x_mean)
        num += x_res * (y[i] - y_mean)
        nom += x_res * x_res
    end
    slope::T = num / nom
    offset::T = y_mean - slope * x_mean
    return offset, slope
end

function baseline_correction1(dir::AbstractString)
    files = readdir(dir)
    for i=eachindex(files)
        lhd = LHDataStore(joinpath(dir, files[i]))
        samples = lhd["ICPC/samples"][:]
        stats = signalstats.(samples, 1, 1800)
        samples = shift_waveform.(samples, -stats.mean)
    end
end