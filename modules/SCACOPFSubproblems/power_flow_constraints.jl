# function to add a rectangular power flow constraint to a model

function addpowerflowcon!(m::JuMP.Model, pq::Union{JuMP.VariableRef, Real},
                          vi::Union{JuMP.VariableRef, Real},
                          vj::Union{JuMP.VariableRef, Real},
                          thetai::Union{JuMP.VariableRef, Real},
                          thetaj::Union{JuMP.VariableRef, Real},
                          A::Real, B::Real, C::Real, Theta::Real=0)
    return @NLconstraint(m, pq == A*vi^2 + B*vi*vj*cos(thetai - thetaj + Theta) + 
                         C*vi*vj*sin(thetai - thetaj + Theta))
end

## function to add all power flow constraints to a model

function addpowerflowcons!(m::JuMP.Model,
                           v_n::AbstractVector, theta_n::AbstractVector,
                           p_li::AbstractArray, q_li::AbstractArray,
                           p_ti::AbstractArray, q_ti::AbstractArray,
                           b_s::AbstractVector,
                           p_g::AbstractVector, q_g::AbstractVector,
                           pslackm_n::AbstractVector, pslackp_n::AbstractVector,
                           qslackm_n::AbstractVector, qslackp_n::AbstractVector,
                           sslack_li::AbstractArray, sslack_ti::AbstractArray,
                           psd::SCACOPFdata, con::GenericContingency=GenericContingency();
                           aux_slack_i::Union{Nothing, AbstractArray}=nothing)
    
    # check whether we are in base case or a contingency
    if isequal_struct(con, GenericContingency())
        RateSymb = :RateBase
    else
        RateSymb = :RateEmer
    end
    
    # thermal limits and power flows -- lines
    for l=1:nrow(psd.L)
        
        # check if line is out
        if length(con.lines_out) > 0 && l in con.lines_out
            JuMP.fix.(p_li[l,:], 0.0, force=true)
            JuMP.fix.(q_li[l,:], 0.0, force=true)
            JuMP.fix.(sslack_li[l,:], 0.0, force=true)
            continue
        end
        
        # check whether branch is short circuit
        y_l = psd.L[l, :G] + im * psd.L[l, :B]
        if isinf(abs(y_l))
            # short circuit branch constraints
            @constraint(m, p_li[l, 1] + p_li[l, 2] == 0)
            @constraint(m, q_li[l, 1] + q_li[l, 2] == 0)
            from = psd.L_Nidx[l,1]
            to = psd.L_Nidx[l,2]
            @constraint(m, v_n[from] == v_n[to])
            @constraint(m, theta_n[from] == theta_n[to])
        else
            # non linear power flow constraints
            for i = 1:2
                from = psd.L_Nidx[l,i]
                to = psd.L_Nidx[l,3-i]
                addpowerflowcon!(m, p_li[l,i],
                                 v_n[from], v_n[to], theta_n[from], theta_n[to],
                                 psd.L[l,:G], -psd.L[l,:G], -psd.L[l,:B])
                addpowerflowcon!(m, q_li[l,i],
                                 v_n[from], v_n[to], theta_n[from], theta_n[to],
                                 -psd.L[l,:B] - psd.L[l,:Bch]/2, psd.L[l,:B],
                                 -psd.L[l,:G])
            end
        end
        
        # thermal limit constraints
        if !isinf(psd.L[l,RateSymb])
            for i = 1:2
                from = psd.L_Nidx[l,i]
                if isnothing(aux_slack_i)
                    @constraint(m, p_li[l,i]^2 + q_li[l,i]^2 <=
                                (psd.L[l,RateSymb]*v_n[from] + sslack_li[l,i])^2)
                else
                    eq_slack = @variable(m, lower_bound = 0)
                    relax_slack = @variable(m); push!(aux_slack_i, relax_slack);
                    @constraint(m, p_li[l,i]^2 + q_li[l,i]^2 + eqslack -
                                (psd.L[l,RateSymb]*v_n[from] + sslack_li[l,i])^2 ==
                                relax_slack)
                end
            end
        end
        
    end
    
    # thermal limits and power flows -- transformers
    for t=1:size(psd.T, 1)
        
        # check if transformer is out
        if length(con.transformers_out) > 0 && t in con.transformers_out
            JuMP.fix.(p_ti[t,:], 0.0, force=true)
            JuMP.fix.(q_ti[t,:], 0.0, force=true)
            JuMP.fix.(sslack_ti[t,:], 0.0, force=true)
            continue
        end
        
        # define to and from buses
        from = psd.T_Nidx[t,1] 
        to = psd.T_Nidx[t,2] 
        
        # check whether transformer is short circuit
        y_t = psd.T[t, :G] + im * psd.T[t, :B]
        if isinf(abs(y_t))
            # short circuit transformer constraints
            @constraint(m, p_ti[t, 1] + p_ti[t, 2] == 0)
            @constraint(m, q_ti[t, 1] + q_ti[t, 2] == 0)
            @constraint(m, v_n[from]/psd.T[t,:Tau] == v_n[to])
            @constraint(m, theta_n[from] - psd.T[t,:Theta] == theta_n[to])
        else
            # non linear power flow constraints
            addpowerflowcon!(m, p_ti[t,1],
                             v_n[from], v_n[to], theta_n[from], theta_n[to],
                             psd.T[t,:G]/psd.T[t,:Tau]^2 + psd.T[t,:Gm],
                             -psd.T[t,:G]/psd.T[t,:Tau],
                             -psd.T[t,:B]/psd.T[t,:Tau],
                             -psd.T[t,:Theta])
            addpowerflowcon!(m, q_ti[t,1],
                             v_n[from], v_n[to], theta_n[from], theta_n[to],
                             -psd.T[t,:B]/psd.T[t,:Tau]^2 - psd.T[t,:Bm],
                             psd.T[t,:B]/psd.T[t,:Tau],
                             -psd.T[t,:G]/psd.T[t,:Tau],
                             -psd.T[t,:Theta])
            addpowerflowcon!(m, p_ti[t,2],
                             v_n[to], v_n[from], theta_n[to], theta_n[from],
                             psd.T[t,:G], -psd.T[t,:G]/psd.T[t,:Tau],
                             -psd.T[t,:B]/psd.T[t,:Tau], psd.T[t,:Theta])
            addpowerflowcon!(m, q_ti[t,2],
                             v_n[to], v_n[from], theta_n[to], theta_n[from],
                             -psd.T[t,:B], psd.T[t,:B]/psd.T[t,:Tau],
                             -psd.T[t,:G]/psd.T[t,:Tau], psd.T[t,:Theta])
        end
        
        # thermal limit constraints
        if !isinf(psd.T[t, RateSymb])
            if isnothing(aux_slack_i)
                @constraint(m, [i=1:2], p_ti[t,i]^2 + q_ti[t,i]^2 <=
                            (psd.T[t,RateSymb] + sslack_ti[t,i])^2)
            else
                for i = 1:2
                    eq_slack = @variable(m, lower_bound = 0)
                    relax_slack = @variable(m); push!(aux_slack_i, relax_slack);
                    @constraint(m, [i=1:2], p_ti[t,i]^2 + q_ti[t,i]^2 + eq_slack -
                                (psd.T[t,RateSymb] + sslack_ti[t,i])^2 == relax_slack)
                end
            end
        end
        
    end
    
    # balance
    for n = 1:size(psd.N, 1)
        Lidxn = psd.Lidxn[n]
        Lin = psd.Lin[n]
        Tidxn = psd.Tidxn[n]
        Tin = psd.Tin[n]
        if isnothing(aux_slack_i)
            p_relax_slack = 0.0
            q_relax_slack = 0.0
        else
            p_relax_slack = @variable(m); push!(aux_slack_i, p_relax_slack);
            q_relax_slack = @variable(m); push!(aux_slack_i, q_relax_slack);
        end
        @constraint(m, sum(p_g[g] for g in psd.Gn[n]) - 
                    psd.N[n,:Pd] - psd.N[n,:Gsh]*v_n[n]^2 -
                    sum(p_li[Lidxn[lix], Lin[lix]] for lix=1:length(Lidxn)) -
                    sum(p_ti[Tidxn[tix], Tin[tix]] for tix=1:length(Tidxn)) ==
                    pslackp_n[n] - pslackm_n[n] + p_relax_slack)
        @NLconstraint(m, sum(q_g[g] for g in psd.Gn[n]) - psd.N[n,:Qd] -
                      (-psd.N[n,:Bsh] - sum(b_s[s] for s=psd.SShn[n]))*v_n[n]^2 -
                      sum(q_li[Lidxn[lix], Lin[lix]] for lix=1:length(Lidxn)) -
                      sum(q_ti[Tidxn[tix], Tin[tix]] for tix=1:length(Tidxn)) ==
                      qslackp_n[n] - qslackm_n[n] + q_relax_slack)
    end
end

function addpowerflowcons!(v_n::AbstractVector, theta_n::AbstractVector,
                           p_li::AbstractArray, q_li::AbstractArray,
                           p_ti::AbstractArray, q_ti::AbstractArray,
                           b_s::AbstractVector,
                           p_g::AbstractVector, q_g::AbstractVector,
                           pslackm_n::AbstractVector, pslackp_n::AbstractVector,
                           qslackm_n::AbstractVector, qslackp_n::AbstractVector,
                           sslack_li::AbstractArray, sslack_ti::AbstractArray,
                           psd::SCACOPFdata, con::GenericContingency=GenericContingency();
                           aux_slack_i::Union{Nothing, AbstractArray}=nothing)
    
    # check whether we are in base case or a contingency
    if isequal_struct(con, GenericContingency())
        RateSymb = :RateBase
    else
        RateSymb = :RateEmer
    end
        
    sum_pg = Any[]
    pd = Vector{Float64}(undef, size(psd.N, 1))
    gsh = Vector{Float64}(undef, size(psd.N, 1))
    vn = Vector{Float64}(undef, size(psd.N, 1))
    gshvn2 = Vector{Float64}(undef, size(psd.N, 1))
    sum_p_li = Any[]
    sum_p_ti = Any[]
    p_relax = Vector{Float64}(undef, size(psd.N, 1))
    p_mn = Vector{Float64}(undef, size(psd.N, 1))
    p_pn = Vector{Float64}(undef, size(psd.N, 1))

    sum_qg = Any[]
    qd = Vector{Float64}(undef, size(psd.N, 1))
    Bsh = Vector{Float64}(undef, size(psd.N, 1))
    ssh = Any[]
    sshvn2 = Any[]
    sum_q_li = Any[]
    sum_q_ti = Any[]
    q_relax = Vector{Float64}(undef, size(psd.N, 1))
    q_mn = Vector{Float64}(undef, size(psd.N, 1))
    q_pn = Vector{Float64}(undef, size(psd.N, 1))

    # balance
    for n = 1:size(psd.N, 1)
        Lidxn = psd.Lidxn[n]
        Lin = psd.Lin[n]
        Tidxn = psd.Tidxn[n]
        Tin = psd.Tin[n]
        if isnothing(aux_slack_i)
            p_relax_slack = 0.0
            q_relax_slack = 0.0
        else
            p_relax_slack = @variable(m); push!(aux_slack_i, p_relax_slack);
            q_relax_slack = @variable(m); push!(aux_slack_i, q_relax_slack);
        end

        if length(psd.Gn[n]) == 0
            push!(sum_pg, NaN)
            push!(sum_qg, NaN)
        else
            push!(sum_pg, JuMP.value(sum(p_g[g] for g in psd.Gn[n])))
            push!(sum_qg, JuMP.value(sum(q_g[g] for g in psd.Gn[n])))
        end
        if length(Lidxn) == 0
            push!(sum_p_li, NaN)
            push!(sum_q_li, NaN)
        else
            push!(sum_p_li, JuMP.value(sum(p_li[Lidxn[lix], Lin[lix]] for lix=1:length(Lidxn))))
            push!(sum_q_li, JuMP.value(sum(q_li[Lidxn[lix], Lin[lix]] for lix=1:length(Lidxn)) ))
        end
        if length(Tidxn) == 0
            push!(sum_p_ti, NaN)
            push!(sum_q_ti, NaN)
        else
            push!(sum_p_ti, -JuMP.value(sum(p_ti[Tidxn[tix], Tin[tix]] for tix=1:length(Tidxn))))
            push!(sum_q_ti, -JuMP.value(sum(q_ti[Tidxn[tix], Tin[tix]] for tix=1:length(Tidxn)) ))
        end
        if length(psd.SShn[n]) == 0
            push!(sshvn2, NaN)
            push!(ssh, NaN)
        else
            push!(sshvn2, -JuMP.value((-psd.N[n,:Bsh] - sum(b_s[s] for s=psd.SShn[n]))*v_n[n]^2))
            push!(ssh,  JuMP.value(sum(b_s[s] for s=psd.SShn[n])))
        end

        pd[n] = -psd.N[n,:Pd]
        gsh[n] = psd.N[n,:Gsh]
        vn[n] = JuMP.value(v_n[n])
        gshvn2[n] = -JuMP.value(psd.N[n,:Gsh]*v_n[n]^2)
        p_relax[n] = p_relax_slack
        p_mn[n] = JuMP.value(pslackm_n[n])
        p_pn[n] = JuMP.value(pslackp_n[n])
    
        qd[n] = -psd.N[n,:Qd]
        Bsh[n] = psd.N[n,:Bsh]
        q_relax[n] = q_relax_slack
        q_mn[n] = JuMP.value(qslackm_n[n])
        q_pn[n] = JuMP.value(qslackp_n[n])
    end
    write_pf_constraint_info("./PF_constraint", "/power_constraints", psd, sum_pg, pd, gsh, vn, gshvn2, sum_p_li,
                            sum_p_ti, p_relax, sum_qg, qd, Bsh, ssh, sshvn2, sum_q_li, sum_q_ti, q_relax,
                            p_mn, p_pn, q_mn, q_pn)
end

function write_pf_constraint_info(OutDir::String, filename::String, psd::SCACOPFdata, 
                                sum_pg::Array{Any}, pd::Array{Float64}, gsh::Array{Float64}, 
                                v_n::Array{Float64},  gshvn2::Array{Float64}, sum_p_li::Array{Any},
                                sum_p_ti::Array{Any},  p_relax::Array{Float64}, 
                                sum_qg::Array{Any}, qd::Array{Float64},  Bsh::Array{Float64}, 
                                ssh::Array{Any},  sshvn2::Array{Any}, 
                                sum_q_li::Array{Any}, sum_q_ti::Array{Any}, q_relax::Array{Float64},
                                p_mn::Array{Float64}, p_pn::Array{Float64},
                                q_mn::Array{Float64}, q_pn::Array{Float64})
    display(filename)
    if !ispath(OutDir)
		mkpath(OutDir)
	end
    f = open(OutDir * filename, "w")
    @printf(f, "--Active power constraint\n")
    @printf(f, "Bus, sum_pg, pd, gsh, v_n, gshvn2, sum_p_li, sum_p_ti, pslackm_n, pslackp_n, p_relax\n")
	for n=1:size(psd.N, 1)
		@printf(f, "%d, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n", 
                psd.N[n,:Bus], sum_pg[n], pd[n], gsh[n], v_n[n], gshvn2[n], sum_p_li[n], sum_p_ti[n], 
                p_mn[n], p_pn[n], p_relax[n])
	end
	
	# write line section
    @printf(f, "--Reactive power constraint\n")
    @printf(f, "Bus, sum_qg, qd, Bsh, Ssh, v_n, BshSshvn2, sum_q_li, sum_q_ti, qslackm_n, qslackp_n, q_relax\n")
	for n=1:size(psd.N, 1)
		@printf(f, "%d, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n", 
                psd.N[n,:Bus], sum_qg[n], qd[n], Bsh[n], ssh[n], v_n[n], sshvn2[n], sum_q_li[n], sum_q_ti[n],
                q_mn[n], q_pn[n], q_relax[n])
	end
end
