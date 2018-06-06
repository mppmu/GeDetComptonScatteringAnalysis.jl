# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).


function triangular_dither(x::Real, width::Real = one(typeof(x)))
    T = float(typeof(x))
    r = rand(T)
    tr = (r >= T(0.5)) ? - sqrt(2 - 2*r) + 1 : sqrt(2*r) - 1
    x + tr * width
end


pos_in_mm(raw_x::Real) = raw_x * 1E-6 * 1E-3

edep_in_keV(raw_edep::Real) = raw_edep * 1E-3


function polaris_events2df(events::PolarisEvents)
    T = Float32
    inv_1000 = 1 / T(1000)

    DataFrames.DataFrame(
        evt_no = events.evt_no,
        evt_nhits = events.evt_nhits,
        evt_t = events.evt_t * 1E-9,
        evt_issync = events.evt_issync,
        hit_detno = events.hit_detno,
        hit_x = deepmap(x -> triangular_dither(T(x)) * inv_1000, events.hit_x),
        hit_y = deepmap(x -> triangular_dither(T(x)) * inv_1000, events.hit_y),
        hit_z = deepmap(x -> triangular_dither(T(x)) * inv_1000, events.hit_z),
        hit_edep = deepmap(x -> triangular_dither(T(x)) * inv_1000, events.hit_edep),
        hit_t = deepmap(x -> x * 1E-9, events.hit_t),
    )
end

export polaris_events2df


function correct_timestamps!(dfs::DataFrame...)
    t_sync = map(df -> df[find(identity, df[:evt_issync]), :evt_t], dfs)
    n_sync_events = min(map(length, t_sync)...)
    resize!.(t_sync, n_sync_events)

    all(x -> length(t_sync) - n_sync_events <= 2, t_sync) || error("Number of sync events doesn't match")

    t_offs = dfs[1][1, :evt_t]

    for i in 2:length(dfs)
        # Determine timestamp mapping relative to reference dataframe:
        t_sync_fit = curve_fit((x, p) -> p[1] + p[2] * x, t_sync[i], t_sync[1], [1.0, 1.0])
        sync_p1, sync_p2 = t_sync_fit.param[1], t_sync_fit.param[2]

        # Correct timestamps:
        dfs[i][:evt_t] .= sync_p1 .+ sync_p2 .* dfs[i][:evt_t] .- t_offs
    end
    dfs[1][:evt_t] .= dfs[1][:evt_t] .- t_offs
    dfs
end

export correct_timestamps!


function find_common_events(df::Tuple{DataFrame,DataFrame}, delta_t::Real)
    t = (df[1][:evt_t], df[2][:evt_t])
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

export find_common_events


const m_e = 510.9989

compton_E_out(E_in::Real, θ::Real) = E_in / (1 + E_in/m_e * (1 - cos(θ)))

function compton_theta(E_in::Real, E_out::Real)
    cos_theta = 1 - m_e/E_out + m_e/E_in
    T = typeof(cos_theta)
    (-1 < cos_theta < 1) ? T(acos(cos_theta)) : T(NaN)
end

#=
Tests:
compton_theta(622, 600) * 180 / pi
compton_theta(622, 181.2) * 180 / pi
=#
