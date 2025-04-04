# function to compute primal starting point

function get_primal_starting_point(psd::SCACOPFdata,
                                   sol::Union{Nothing, BasecaseSolution}=nothing)
    x0 = Dict{Symbol, Array{Float64}}()
    if isnothing(sol)
        full_solution = get_full_initial_solution(psd)
        x0[:v_n] = full_solution[1]
        x0[:theta_n] = full_solution[2]
        x0[:b_s] = full_solution[3]
        x0[:p_g] = full_solution[4]
        x0[:q_g] = full_solution[5]
        x0[:p_li] = full_solution[6]
        x0[:q_li] = full_solution[7]
        x0[:p_ti] = full_solution[8]
        x0[:q_ti] = full_solution[9]
        x0[:pslackm_n] = full_solution[10]
        x0[:pslackp_n] = full_solution[11]
        x0[:qslackm_n] = full_solution[12]
        x0[:qslackp_n] = full_solution[13]
        x0[:sslack_li] = full_solution[14]
        x0[:sslack_ti] = full_solution[15]
        x0[:c_g] = full_solution[16]
    else
        full_solution = get_full_solution(psd, sol)
        x0[:v_n] =  sol.v_n
        x0[:theta_n] = sol.theta_n
        x0[:b_s] = sol.b_s
        x0[:p_g] = sol.p_g
        x0[:q_g] = sol.q_g
        x0[:p_li] = full_solution[1]
        x0[:q_li] = full_solution[2]
        x0[:p_ti] = full_solution[3]
        x0[:q_ti] = full_solution[4]
        x0[:pslackm_n] = full_solution[5]
        x0[:pslackp_n] = full_solution[6]
        x0[:qslackm_n] = full_solution[7]
        x0[:qslackp_n] = full_solution[8]
        x0[:sslack_li] = full_solution[9]
        x0[:sslack_ti] = full_solution[10]
        x0[:c_g] = full_solution[11]
    end
    return x0
end

function get_primal_starting_point(psd::SCACOPFdata, con::GenericContingency,
                                   base_sol::BasecaseSolution,
                                   prev_sol::Union{Nothing,
                                                   ContingencySolution})
    x0 = Dict{Symbol, Array{Float64}}()
    if isnothing(prev_sol) || !isequal_struct(prev_sol.cont_alt, con)
        full_solution = get_full_initial_solution(psd, con, base_sol)
        x0[:v_nk] = full_solution[1]
        x0[:theta_nk] = full_solution[2]
        x0[:b_sk] = full_solution[3]
        x0[:p_gk] = full_solution[4]
        x0[:q_gk] = full_solution[5]
        x0[:p_lik] = full_solution[6]
        x0[:q_lik] = full_solution[7]
        x0[:p_tik] = full_solution[8]
        x0[:q_tik] = full_solution[9]
        x0[:pslackm_nk] = full_solution[10]
        x0[:pslackp_nk] = full_solution[11]
        x0[:qslackm_nk] = full_solution[12]
        x0[:qslackp_nk] = full_solution[13]
        x0[:sslack_lik] = full_solution[14]
        x0[:sslack_tik] = full_solution[15]
    else
        @assert isequal_struct(prev_sol.cont_alt, con)
        full_solution = get_full_solution(psd, prev_sol)
        x0[:v_nk] =  prev_sol.v_n
        x0[:theta_nk] = prev_sol.theta_n
        x0[:b_sk] = prev_sol.b_s
        x0[:p_gk] = prev_sol.p_g
        x0[:q_gk] = prev_sol.q_g
        x0[:p_lik] = full_solution[1]
        x0[:q_lik] = full_solution[2]
        x0[:p_tik] = full_solution[3]
        x0[:q_tik] = full_solution[4]
        x0[:pslackm_nk] = full_solution[5]
        x0[:pslackp_nk] = full_solution[6]
        x0[:qslackm_nk] = full_solution[7]
        x0[:qslackp_nk] = full_solution[8]
        x0[:sslack_lik] = full_solution[9]
        x0[:sslack_tik] = full_solution[10]
    end
    return x0
end

function get_primal_starting_point(psd::SCACOPFdata, con::GenericContingency)
    x0 = Dict{Symbol, Array{Float64}}()
    full_solution = get_full_initial_solution(psd, con)
    x0[:v_nk] = full_solution[1]
    x0[:theta_nk] = full_solution[2]
    x0[:b_sk] = full_solution[3]
    x0[:p_gk] = full_solution[4]
    x0[:q_gk] = full_solution[5]
    x0[:p_lik] = full_solution[6]
    x0[:q_lik] = full_solution[7]
    x0[:p_tik] = full_solution[8]
    x0[:q_tik] = full_solution[9]
    x0[:pslackm_nk] = full_solution[10]
    x0[:pslackp_nk] = full_solution[11]
    x0[:qslackm_nk] = full_solution[12]
    x0[:qslackp_nk] = full_solution[13]
    x0[:sslack_lik] = full_solution[14]
    x0[:sslack_tik] = full_solution[15]
    return x0
end