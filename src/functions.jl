# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).
using RadiationDetectorDSP: SamplesOrWaveform, RealQuantity

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

@fastmath function linear_regression(x::Vector{<:Real}, y::Vector{<:Real}
)::Tuple
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

@fastmath function linear_regression(x::AbstractVector{T}, y::AbstractVector{T}
)::Tuple{T, T} where {T <: AbstractFloat}
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

"""
    linreg(y::AbstractVector{T}, y₀::Float64
        )::Tuple{T, T} where {T<:AbstractFloat}

perform linear regression on `log.(max.(y .- y₀, 1))` assuming that 
`x=0:L-1` with `L=length(y)` and `y` being the data.
"""
@fastmath function loglinreg(y::AbstractVector{T}, y₀::Float64) where {T<:Real}
    L = length(y)
    x̄::Float64 = (L-1) / 2.
    ȳ::Float64 = 0.
    @inbounds @simd for i=Base.OneTo(L)
        ȳ += log(max((y[i] - y₀), 1))
    end
    ȳ /= L
    num::Float64 = 0.
    nom::Float64 = 0.
    @inbounds for i=Base.OneTo(L)
        x_res = i - 1 - x̄
        num += x_res * (log(max((y[i] - y₀), 1)) - ȳ)
        nom += x_res * x_res
    end
    slope::Float64 = num / nom
    offset::Float64 = ȳ - slope * x̄
    return offset, slope
end

"""
    get_tau(input::AbstractSamples{T}, baseline_length::Int, step::T, 
    tail_start::Int, tail_stop::Int)::Float64 where {T<:Real}

Determine the decay time tau from the tail and the baseline_mean from 
the beginning of the pulse `input`, where `input` is from type 
`AbstractSamples{T}`. `baseline_length` is the length of the baseline in 
number of samples. `tail_start` and `tail_stop` define the beginning and 
end of the tail, both in number of samples.
"""
function get_tau(y::AbstractVector{T}, baseline_length::Int,
tail_start::Int, tail_stop::Int) where {T<:Real}  
    # @assert tail_stop >= tail_start
    y₀::Float64 = ustrip(mean(@view y[1:baseline_length]))
    fitresult = loglinreg((@view y[tail_start:tail_stop]), y₀)
    (tau = -1 / fitresult[2], baseline_mean=y₀)
end

"""
    get_tau(input::RDWaveform{T}, baseline_length::Int, tail_start::Int, 
    tail_stop::Int)::T where {T<:RealQuantity})

Determine the decay time tau from the tail and the baseline_mean from 
the beginning of the pulse 'input', where `input` is an `RDWaveform`.
"""
function get_tau(x::RDWaveform{T}, baseline_length::RealQuantity, 
tail_start::RealQuantity, tail_stop::RealQuantity) where {T<:RealQuantity} 
    _baseline_length = _get_index(convert(T, baseline_length), x.time)
    _tail_start = _get_index(convert(T, tail_start), x.time)
    _tail_stop = _get_index(convert(T, tail_stop), x.time)
    _tau, baseline_mean = 
        get_tau(x.signal, _baseline_length, _tail_start, _tail_stop)
    (tau = _tau * step(x.time), baseline_mean=T(baseline_mean))
end

@inline cauchy_model(x, p)  = 
    p[3]/pi .* p[1]./(p[1].^2 .+ (x .- p[2]).^2)

@inline function _get_index(stop::U, x_axis::AbstractVector{T}
)::Int where {T<:RealQuantity, U<:RealQuantity}
    first_x, step_x =  first(x_axis), step(x_axis)
    round(Int, ustrip(NoUnits, (stop - first_x) / step_x)) + firstindex(x_axis)
end

@inline findmax(x::AbstractSamples{T}) where {T<:Real} = maximum(x)
@inline findmax(x::RDWaveform{T, U}) where {T<:Real, U<:RealQuantity} = 
    findmax(x.signal)