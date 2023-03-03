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
        Atmosphere, CreateAtmosphere,
        Surface, CreateSurface,
        SimpleIslandModel, CreateModel,
        calculateFₐ, calculateFₛ,
        run, solveanalytic

"""
    Atmosphere

Abstract supertype for the different atmospheric layer types
"""
abstract type Atmosphere end

"""
    Surface

Abstract supertype for the different surface layer types
"""
abstract type Surface end

struct Variables{FT<:Real}
    t :: Vector{FT}
   S₀ :: Vector{FT}
   Tₛ :: Vector{FT}
   Tₐ :: Vector{FT}
   Fₛ :: Vector{FT}
   Fₐ :: Vector{FT}
   temp :: Vector{FT} # temp = [t,S₀,Tₛ,Tₐ,Fₛ,Fₐ,δt]
   stat :: Vector{FT} # temp = [t,S₀,Tₛ,Tₐ,Fₛ,Fₐ]
end

const σ = 5.670374419e-8

modulelog() = "$(now()) - SimpleIslandModels.jl"

include("model/surface.jl")
include("model/atmosphere.jl")
include("model/model.jl")

include("variables.jl")

include("calculate/basics.jl")
include("calculate/atmosphere.jl")
include("calculate/surface.jl")
include("calculate/cloud.jl")
include("calculate/analytic.jl")

include("run.jl")
include("analytic.jl")

end # module SimpleIslandModels
