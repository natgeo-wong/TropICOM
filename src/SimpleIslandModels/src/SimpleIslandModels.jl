module SimpleIslandModels

## Base Modules Used
using Dates
using DelimitedFiles
using Logging
using Printf
using Statistics

## Modules Used
using ProgressMeter

import Base: eltype, show, run

export
        CreateFixedSSTModel, CreateFixedAtmosModel, CreateMixedLayerModel

"""
    SimpleIslandModel

Abstract supertype for different models that YangShallowWater can run.
"""
abstract type SimpleIslandModel end

struct Variables{FT<:Real}
    t :: Vector{FT}
   S₀ :: Vector{FT}
   Tₛ :: Vector{FT}
   Tₐ :: Vector{FT}
   temp :: Vector{FT} # temp = [t,S₀,Tₛ,Tₐ,δt]
   stat :: Vector{FT} # temp = [t,S₀,Tₛ,Tₐ,δt]
end

const σ = 5.670374419e-8

modulelog() = "$(now()) - SimpleIslandModels.jl"

include("FixedAtmos.jl")
include("FixedSST.jl")
include("MixedLayer.jl")

include("variables.jl")
include("calculate.jl")
include("run.jl")

end # module SimpleIslandModels
