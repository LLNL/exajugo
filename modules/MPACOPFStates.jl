# this is part of module MPACOPF

mutable struct MPBaseState
    v_n::JuMP.Containers.DenseAxisArray
    theta_n::JuMP.Containers.DenseAxisArray
    b_s::JuMP.Containers.DenseAxisArray
    p_g::JuMP.Containers.DenseAxisArray
    q_g::JuMP.Containers.DenseAxisArray
end

mutable struct biddingRes
    prodcost::Float64
    penalty_base::Float64
    plantrev::Float64
    # plantpower::JuMP.Containers.DenseAxisArray
    plantpower::Array{Float64}
    # LMPs::Array{Float64}
    LMPs::JuMP.Containers.DenseAxisArray{Float64}
    LMPsN::JuMP.Containers.DenseAxisArray{Float64}
end

mutable struct TranLoss
    LineLoss::Union{NTuple, Vector{Float64}}
    GenTotal::Union{NTuple, Vector{Float64}}
    Eta::Union{NTuple, Vector{Float64}}
    p_line::JuMP.Containers.DenseAxisArray
    p_tran::JuMP.Containers.DenseAxisArray
    Gentotbus::JuMP.Containers.DenseAxisArray
    Llossbus::JuMP.Containers.DenseAxisArray  
    GenCostStep::Vector{Float64}
    Lactive::DataFrame
    Tactive::DataFrame
end

mutable struct slackv
    spm::Array{Float64}
    spp::Array{Float64}
    sqm::Array{Float64}
    sqp::Array{Float64}
    sl::Array{Float64,3}
    st::Array{Float64,3}

end
mutable struct MPACOPFBaseState
    #base
    base_results::MPBaseState

    #bidding solutions
    bidding_results::biddingRes

    #slack variables
    slackn::slackv

    #transmission losses
    TLoss::TranLoss
end    

mutable struct  Ylinkdecomts
    YUk::Union{NTuple, Vector{Float64}}
    YZk::Union{NTuple, Vector{Float64}}
    YVk::Union{NTuple, Vector{Float64}}
    YUkm::Union{NTuple, Vector{Float64}}
    YVkm::Union{NTuple, Vector{Float64}}
    ts::Int64   

end

function MPACOPFState(base::MPBaseState, bidding::biddingRes, slackness::slackv, Loss::TranLoss)

return MPACOPFBaseState(base, bidding, slackness, Loss)
end
