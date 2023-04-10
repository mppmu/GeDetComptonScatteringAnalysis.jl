# This file is a part of GeDetComptonScatteringAnalysis.jl, licensed under the MIT License (MIT).

using Unitful: 𝐋, 𝐓, Quantity, lookup_units

const eV = typeof(lookup_units([Unitful], :eV))
const ns = typeof(lookup_units([Unitful], :ns))
const μm = typeof(lookup_units([Unitful], :μm))
const mm = typeof(lookup_units([Unitful], :mm))

const QuantityMM{T} = Quantity{T, 𝐋, mm}
const MaybeWithUnitsMM{T} = Union{T, QuantityMM{T}}


const cztTuple = NamedTuple{
    (:evt_no, :evt_t, :evt_nhits, :evt_issync, :hit_edep, :hit_t, :hit_detno, :hit_x, :hit_y, :hit_z), 
    Tuple{
        Int32, # evt_no
        Quantity{Int64, 𝐓, ns}, # evt_t
        Int32, # evt_nhits
        Bool, # evt_issync
        Vector{Quantity{Int32, dimension(1u"eV"), eV}}, # hit_edep
        Vector{Quantity{Int64, 𝐓, ns}}, # hit_t
        Vector{Int32}, # hit_detno
        Vector{Quantity{Int32, 𝐋, μm}}, # hit_x
        Vector{Quantity{Int32, 𝐋, μm}}, # hit_y
        Vector{Quantity{Int32, 𝐋, μm}}  # hit_z
    }
}

const cztTable = TypedTables.Table{
    cztTuple, 1,
    NamedTuple{
        (:evt_no, :evt_t, :evt_nhits, :evt_issync, :hit_edep, :hit_t, :hit_detno, :hit_x, :hit_y, :hit_z), 
        Tuple{
            Vector{Int32}, # evt_no
            Vector{Quantity{Int64, 𝐓, ns}}, # evt_t
            Vector{Int32}, # evt_nhits
            BitVector, # evt_issync
            VectorOfVectors{Quantity{Int32, dimension(1u"eV"), eV}, Vector{Quantity{Int32, dimension(1u"eV"), eV}}, Vector{Int64}, Vector{Tuple{}}}, # hit_edep
            VectorOfVectors{Quantity{Int64, 𝐓, ns}, Vector{Quantity{Int64, 𝐓, ns}}, Vector{Int64}, Vector{Tuple{}}}, # hit_t
            VectorOfVectors{Int32, Vector{Int32}, Vector{Int64}, Vector{Tuple{}}}, # hit_detno
            VectorOfVectors{Quantity{Int32, 𝐋, μm}, Vector{Quantity{Int32, 𝐋, μm}}, Vector{Int64}, Vector{Tuple{}}}, # hit_x
            VectorOfVectors{Quantity{Int32, 𝐋, μm}, Vector{Quantity{Int32, 𝐋, μm}}, Vector{Int64}, Vector{Tuple{}}}, # hit_y
            VectorOfVectors{Quantity{Int32, 𝐋, μm}, Vector{Quantity{Int32, 𝐋, μm}}, Vector{Int64}, Vector{Tuple{}}}  # hit_z
        }
    }
}

const detTuple = NamedTuple{
    (:evt_no, :chid, :evt_t, :DAQ_energy, :samples), 
    Tuple{
        Int32, # evt_no
        Int32, # chid
        Quantity{Int64, 𝐓, ns}, # evt_t
        Int32, # DAQ_energy
        Vector{Int16} # samples
    }
}

const detTable = TypedTables.Table{
    detTuple, 1, 
    NamedTuple{
        (:evt_no, :chid, :evt_t, :DAQ_energy, :samples), 
        Tuple{
            Vector{Int32}, # evt_no
            Vector{Int32}, # chid
            Vector{Quantity{Int64, 𝐓, ns}}, # evt_t
            Vector{Int32}, # DAQ_energy
            VectorOfVectors{Int16, Vector{Int16}, Vector{Int64}, Vector{Tuple{}}} # samples
        }
    }
}

const detTable2 = AbstractVector{
    NamedTuple{
        (:evt_no, :chid, :evt_t, :DAQ_energy, :samples), 
        Tuple{
            Int32, 
            Int32, 
            Quantity{Int64, 𝐓, ns}, 
            Int32, 
            Vector{Int16}
        }
    }
}

_detTable() = 
    TypedTables.Table(
        evt_no = Vector{Int32}(), # evt_no
        chid = Vector{Int32}(), # chid
        evt_t = Vector{Quantity{Int64, 𝐓, ns}}(), # evt_t
        DAQ_energy = Vector{Int32}(), # DAQ_energy
        samples = VectorOfVectors{Int16}() # samples
    )
    
_cztTable() = 
    TypedTables.Table(
        evt_no = Vector{Int32}(), # evt_no
        evt_t = Vector{Quantity{Int64, 𝐓, ns}}(), # evt_t
        evt_nhits = Vector{Int32}(), # evt_nhits
        evt_issync = BitVector(), # evt_issync
        hit_edep = VectorOfVectors{Quantity{Int32, dimension(1u"eV"), eV}}(), # hit_edep
        hit_t = VectorOfVectors{Quantity{Int64, 𝐓, ns}}(), # hit_t
        hit_detno = VectorOfVectors{Int32}(), # hit_detno
        hit_x = VectorOfVectors{Quantity{Int32, 𝐋, μm}}(), # hit_x
        hit_y = VectorOfVectors{Quantity{Int32, 𝐋, μm}}(), # hit_y
        hit_z = VectorOfVectors{Quantity{Int32, 𝐋, μm}}(), # hit_z
    )
export _detTable, _cztTable