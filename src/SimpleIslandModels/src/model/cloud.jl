struct CloudAlbedo
    calculateαₐ :: Function
    params :: CloudAlbedoType
end

struct DummyCloud{FT<:Real} <: CloudAlbedoType
    δα :: FT
end

DummyCloud(FT=Float64; δa::Real=0) = DummyCloud{FT}(δa)

dummycalculateαₐ(
    αₐ::Real, Ts::Real, Ta::Real, δt::Real,
    cld::DummyCloud, atm::Atmosphere, sfc::Surface
) = αₐ + cld.δα

"""
To create a new CloudAlbedo parameterization, one must create a CloudAlbedoType
Type, which contains the parameters, a function to input the chosen paramters,
AND a corresponding function to calculate the cloud-albedo response.

For example, if we wanted to calculate a simple cloud-albedo scheme, we can do
```
struct SimpleCloudExample{FT<:Real} <: CloudAlbedoType
    cₛα :: FT
    cₐα :: FT
    mαₐ :: FT
end

SimpleCloudExample(FT = Float64;
    cₛα :: Real = 0.001
    cₐα :: Real = 0.001
    mαₐ :: Real = 0.8
) = SimpleCloudExample{FT}(cₛα,cₐα,mαₐ)

function calculateαₐ_example(
    αₐ::Real, Ts::Real, Ta::Real, δt::Real,
    cld::SimpleCloudExample, atm::Atmosphere, sfc::Surface
)

    dTs = Ts - sfc.Tsr
    dTa = Ta - atm.Tar

    αₐ += (dTs * m.cₛα + dTa * m.cₐα) * δt

    if αₐ<=0; αₐ = 0 end
    if αₐ>=m.mαₐ; αₐ = m.mαₐ end

    return αₐ

end
```

Then, you should be able to do
```
cld = CloudAlbedo(calculateαₐ_example,SimpleCloudExample())
```
"""