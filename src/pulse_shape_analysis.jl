# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

# TODO: replace pulse-shape processing functions by functions from RadiationDetectorDSP.jl

function normalizewf!(x::AbstractSamples{T}; tail_length::Int=300
) where {T <: AbstractFloat}
    tail = mean(@view x[end-tail_length:end])
    x ./= tail
end

function normalizewf!(x::RDWaveform{U, T}; kwargs...) where {U, T}
    normalizewf!(x.signal; kwargs...)
    x
end

function normalizewf(X::AbstractVector{<:AbstractVector{T}}; 
tail_length::Int=300) where {T <: AbstractFloat}
    y = zero(X[1])
    max = 0.
    @inbounds for i=eachindex(X)
        max += mean(@view X[i][end-tail_length:end])
    end
    @inbounds for i=eachindex(X)
        y .+= X[i] ./ max
    end
    y
end

# TODO: check if ! functions would give better performance 
# (certainly better allocs)
# baseline_corr!(x::AbstractVector{T}; m::Int=500) where {T} = begin
#     x̅ = mean(x[1:m])
#     x .-= x̅
# end

baseline_corr!(x::AbstractVector{T}; m::Int=500) where {T} = begin
    x̅ = mean(@view x[1:m])
    x .- x̅
end

# TODO: check if @fastmath really has an advantage
"""
    time_align(wf::AbstractVector{T}; 
        p::T=0.5, window::Tuple{Int, Int}=(500, 500), l::Int=300
    ) where {T <: AbstractFloat}

Find a window of size `sum(window) + 1` for the waveform `wf` such that 
the point where the `p`-th fraction of the maximum value is reached is 
at the center of the window. The maximum value is determined by the last 
`l` values.
"""
function time_align(wf::AbstractVector{T}; p::Float64=0.5, 
window::Tuple{Int, Int}=(500, 500), l::Int=300) where {T<:AbstractFloat}
    @assert l < length(wf) "tail length l cant be longer than the waveform"
    m = mean(@view wf[end-l:end]) * p
    time_align(wf, m; window=window)
end

@fastmath @inbounds function time_align(wf::AbstractSamples{T}, m::T;
window::Tuple{Int, Int}=(500, 500)) where {T<:AbstractFloat}
    @assert length(wf) > sum(window) + 1 "waveform must be longer than windows"
    inter = Intersect(100.)
    idx = floor(Int, inter(wf, m).x)
    (idx ≤ window[1] || length(wf) - window[2] < idx) && return @view wf[:]
    return @view wf[idx-window[1]:idx+window[2]]
end
    
@fastmath function decay_correction(wf::AbstractVector{T}, tau::T
) where {T<:AbstractFloat}
    rp = Vector{T}(undef, length(wf))
    rp[1] = wf[1]
    @inbounds for i in 2:length(rp)
        rp[i] = wf[i] + rp[i-1] - wf[i-1] * tau
    end
    rp
end

# TODO: replace by appropriate filters from DSP/RadiationDetectorDSP
function smooth20(wv::AbstractVector{T}) where {T <: AbstractFloat}
    w::Vector{T} = T[ 0.0,
                     0.010329700377119983,
                     0.022145353094771694,
                     0.034827599769728275,
                     0.04766125969628491,
                     0.05988618988841385,
                     0.07075302574183234,
                     0.07957926681837389,
                     0.08580108206476868,
                     0.0890165225487064,
                     0.0890165225487064,
                     0.08580108206476868,
                     0.07957926681837389,
                     0.07075302574183234,
                     0.05988618988841385,
                     0.04766125969628491,
                     0.034827599769728275,
                     0.022145353094771694,
                     0.010329700377119983,
                     0.0 ]
    DSP.filtfilt(w, T[1], wv)
end

@fastmath function find_crossing(wf::AbstractVector{T}, thold::T; 
window::Tuple{Int, Int} = (100, 1100), interpolate::Bool = true
)::T where {T<:AbstractFloat}
    @inbounds begin
        cross = -1
        i = window[1]
        while cross == -1 && i < window[2]
            if wf[i] ≤ thold && wf[i+1] ≥ thold
                cross = i
            else
                i = i+1
            end
        end
        if cross != -1 && interpolate
            cross = cross + (thold - wf[cross])/(wf[cross+1] - wf[cross])     
        end
        return cross
    end
end

function chi_sq_wfs(x::AbstractVector{T}, y::AbstractVector{T}; 
window=400:600, l=150) where {T}
    varx = var(x[1:l])
    vary = var(y[1:l])
    sum((x[window] - y[window]).^2) / ((length(window) - 1)*(varx + vary))
end