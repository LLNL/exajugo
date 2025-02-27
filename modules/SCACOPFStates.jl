# this is part of module SCACOPF

mutable struct ContingencyState
    cont_id::Int      # global index among all contingencies
    cont_idout::Int   # id of the element out
    cont_type::Symbol

    delta::Float64

    v_n::Array{Float64}
    theta_n::Array{Float64}
    b_s::Array{Float64}
    p_g::Array{Float64}
    q_g::Array{Float64}

    # ACOPF related variables, e.g.,
    # p_li, q_li, p_ti, q_ti, pslackm_n, pslackp_n, qslackm_n, qslackp_n, sslack_li, sslack_ti
    # can/should be added below

    penalty::Float64   # JuMP.value(conpenalty)
    objective::Float64 # JuMP.objective_value(m)
end

mutable struct BaseState
    v_n::Array{Float64}
    theta_n::Array{Float64}
    b_s::Array{Float64}
    p_g::Array{Float64}
    q_g::Array{Float64}

    # ACOPF related variables, e.g.,
    # p_li, q_li, p_ti, q_ti, pslackm_n, pslackp_n, qslackm_n, qslackp_n, sslack_li, sslack_ti
    # can/should be added below

    cost::Float64
    penalty::Float64
end

mutable struct SCACOPFState
    prodcost::Float64
    penalty_base::Float64
    penalty_cont::Float64

    #base
    base_state::BaseState

    # contingency solution: a dictionary that maps between a contingency number and the 
    # corresponding contingency ACOPF solution
    #
    # Note: when running in distributed-memory mode, sol_cont's keys are only the 
    # contingency indexes local to the hosting MPI rank
    cont_states::OrderedDict{Int,ContingencyState}
end    

function SCACOPFState(base::BaseState)
	nullarr=Array{Float64, 2}(undef, 0, 0)
	return SCACOPFState(base.cost, base.penalty, 0.0, base, Int[], Int[], Symbol[],
			    nullarr, nullarr, nullarr, nullarr, nullarr, Float64[])
end

function SCACOPFState(prodcostvalue::Float64, penaltybasevalue::Float64, penaltyconvalue::Float64,
                      base::BaseState,
                      contingencies_id::Array{Int,1},
                      contingencies_idout::Array{Int,1},
                      contingencies_type::Array{Symbol,1},
                      v_nk::Array{Float64,2},
                      theta_nk::Array{Float64,2},
                      b_sk::Array{Float64,2},
                      p_gk::Array{Float64,2},
                      q_gk::Array{Float64,2}, 
                      delta_k::Array{Float64,1})
    contingencies = OrderedDict()
    for k=1:size(contingencies_id,1)
        contingencies[contingencies_id[k]] = 
            ContingencyState(contingencies_id[k], contingencies_idout[k], contingencies_type[k],
                             delta_k[k], v_nk[:,k], theta_nk[:,k],
                             b_sk[:,k], p_gk[:,k], q_gk[:,k],
                             0., 0.)
        
    end
    return SCACOPFState(prodcostvalue, penaltybasevalue, penaltyconvalue, base, contingencies)
end

function num_cont(s::SCACOPFState) 
    return length(s.cont_states)
end

function getBaseStates(s::SCACOPFState)
    return s.base_state.v_n, s.base_state.theta_n, s.base_state.b_s, s.base_state.p_g, s.base_state.q_g
end

function getContingencyStates(s::SCACOPFState, kid::Int64)
    return s.cont_states[kid].delta, s.cont_states[kid].v_n, s.cont_states[kid].theta_n, s.cont_states[kid].b_s, s.cont_states[kid].p_g, s.cont_states[kid].q_g
end

function getBuffersSizes(c::ContingencyState)
    return size(c.v_n,1), size(c.theta_n,1), size(c.b_s,1), size(c.p_g,1), size(c.q_g, 1)
end

# return number of doubles necessary to store c 
#   - ints will be stored as doubles)
#   - 0., 1., 2. -> cont_type::Symbol :Generator, :Line, :Transformer
function getTotalSize(c::ContingencyState)
    nv, ntheta, nbs, n_pg, n_qg = getBuffersSizes(c)
    return  nv + ntheta + nbs + n_pg + n_qg + 6
    #cont_id::Int      # global index among all contingencies
    #cont_idout::Int   # id of the element out
    #cont_type::Symbol
    #delta::Float64
    #penalty
    #objective
end

function toArray(c::ContingencyState)
    n = getTotalSize(c)
    nv, ntheta, nbs, npg, nqg = getBuffersSizes(c)
    ret = Array{Float64}(undef,n);
    ret[1] = c.cont_id
    ret[2] = c.cont_idout
    ret[3] = 0.
    if c.cont_type==:Line; ret[3]=1.; end
    if c.cont_type==:Transformer; ret[3]=2.; end

    ret[4] = c.delta
    
    s = 5; e=4+nv
    ret[s:e] = c.v_n

    s = e+1; e=e+ntheta
    ret[s:e] = c.theta_n

    s = e+1; e = e+nbs
    ret[s:e] = c.b_s

    s = e+1; e = e+npg
    ret[s:e] = c.p_g

    s = e+1; e = e+nqg
    ret[s:e] = c.q_g

    s = e+1; e=e+1
    ret[s] = c.penalty

    s = e+1; e=e+1
    ret[s] = c.objective
    
    @assert s==n
    return ret
end

function fromArray(vals::Array{Float64}, t::ContingencyState)
    n = getTotalSize(t)
    nv, ntheta, nbs, npg, nqg = getBuffersSizes(t)

    con_type=:Generator
    if vals[3]==1.; con_type=:Line; end
    if vals[3]==2.; con_type=:Transformer; end
    
   return  ContingencyState(convert(Int64, vals[1]), convert(Int64, vals[2]), con_type,
                            vals[4], # delta
                            vals[5:4+nv], #v
                            vals[5+nv:4+nv+ntheta],
                            vals[5+nv+ntheta:4+nv+ntheta+nbs],
                            vals[5+nv+ntheta+nbs:4+nv+ntheta+nbs+npg],
                            vals[5+nv+ntheta+nbs+npg:4+nv+ntheta+nbs+npg+nqg],
                            vals[end-1],
                            vals[end])
    
end
    # cont_id::Int      # global index among all contingencies
    # cont_idout::Int   # id of the element out
    # cont_type::Symbol

    # delta::Float64

    # v_n::Array{Float64}
    # theta_n::Array{Float64}
    # b_s::Array{Float64}
    # p_g::Array{Float64}
    # q_g::Array{Float64}

    # # ACOPF related variables, e.g.,
    # # p_li, q_li, p_ti, q_ti, pslackm_n, pslackp_n, qslackm_n, qslackp_n, sslack_li, sslack_ti
    # # can/should be added below

    # penalty::Float64   # JuMP.value(conpenalty)
    # objective::Float64 # JuMP.objective_value(m)
