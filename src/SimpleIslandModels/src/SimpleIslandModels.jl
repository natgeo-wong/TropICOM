module SimpleIslandModels

## Base Modules Used
using Dates
using DelimitedFiles
using Logging
using Printf
using Statistics

## Modules Used
using ProgressLogging

import Base: eltype, show, run

modulelog() = "$(now()) - SimpleIslandModels.jl"

end # module SimpleIslandModels
