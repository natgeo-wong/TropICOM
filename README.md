# **<div align="center">TropICOM</div>**

<p align="center">TROPical Island Convection in Observations and Models</p>

**Authored By:** 
* Nathanael Wong (n.wong@nyu.edu)

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> TropICOM

## Related Repositories:

* [SelfAggregation](https://github.com/natgeo-wong/SelfAggregation)

  This repository contains the spinup prm and lsf files for large-domain RCE
  simulations where convective self-aggregation occurs. While the main focus
  of the repository is exploring the science behind self-aggregation of
  convection, some of the spinups here are used in TropICOM as well.

# Installation

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project.  To (locally) reproduce this project, do the
following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> ] activate .
    Activating environment at `path/to/TropICOM`

   (TropICOM) pkg> instantiate
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box.

## Other Acknowledgements
> Project Repository Template generated using [DrWatson.jl](https://github.com/JuliaDynamics/DrWatson.jl) created by George Datseris.