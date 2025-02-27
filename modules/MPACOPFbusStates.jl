# this is part of module MPACOPFbus

mutable struct MPbusSol
    plantpowerAll::JuMP.Containers.DenseAxisArray{Float64}
    LMPsbusSingle::JuMP.Containers.DenseAxisArray{Float64}
    GenAll::Vector{Float64}
    GenCostStep::Vector{Float64}

end

function MPbusState(PplantAll::JuMP.Containers.DenseAxisArray{Float64},
    LMPbus::JuMP.Containers.DenseAxisArray{Float64},
    GenAllstep::Vector{Float64},
    GenCostAll::Vector{Float64})

    return MPbusSol(PplantAll, LMPbus, GenAllstep, GenCostAll)
end