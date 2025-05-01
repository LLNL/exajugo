# function to write power flow block

function write_power_flow_block(io::IO, psd::SCACOPFdata, 
                                p_li::Matrix{Float64}, p_ti::Matrix{Float64})
	
	# write line section
	@printf(io, "--line section\n")
    @printf(io, "From, To, p_from(MW), p_to(MW)\n")
	for l = 1:nrow(psd.L)
		@printf(io, "%d, %d, %.10f, %.10f\n", psd.L[l,:From],
            psd.L[l,:To], psd.MVAbase*p_li[l, 1], psd.MVAbase*p_li[l, 2])
	end    

    # write transformer section
	@printf(io, "--transformer section\n")
    @printf(io, "From, To, p_from(MW), p_to(MW)\n")
	for t = 1:nrow(psd.T)
		@printf(io, "%d, %d, %.10f, %.10f\n", psd.T[t,:From],
            psd.T[t,:To], psd.MVAbase*p_ti[t, 1], psd.MVAbase*p_ti[t, 2])
	end
    return nothing
    
end

# function to write system slack variable block

function write_slack_block(io::IO, psd::SCACOPFdata, 
            pslackm_n::Array{Float64}, pslackp_n::Array{Float64}, 
            qslackm_n::Array{Float64}, qslackp_n::Array{Float64},
            sslack_li::Matrix{Float64}, sslack_ti::Matrix{Float64})

    # write node section
	@printf(io, "--node slack section\n")
    @printf(io, "Bus, pslackm(MW), pslackp(MW), qslackm(MW), qslackp(MW)\n")
	for n=1:size(psd.N, 1)
		@printf(io, "%d, %.10f, %.10f, %.10f, %.10f\n", psd.N[n,:Bus], 
                psd.MVAbase*pslackm_n[n], psd.MVAbase*pslackp_n[n],
                psd.MVAbase*qslackm_n[n], psd.MVAbase*qslackp_n[n])
	end
	
	# write line section
	@printf(io, "--line slack section\n")
    @printf(io, "From, To, sslack_from(MW), sslack_to(MW)\n")
	for l=1:nrow(psd.L)
		@printf(io, "%d, %d, %.10f, %.10f\n", psd.L[l,:From],
            psd.L[l,:To], psd.MVAbase*sslack_li[l, 1], psd.MVAbase*sslack_li[l, 2])
	end    

	# write transformer section
	@printf(io, "--transformer slack section\n")
    @printf(io, "From, To, sslack_from(MW), sslack_to(MW)\n")
	for t=1:nrow(psd.T)
		@printf(io, "%d, %d, %.10f, %.10f\n", psd.T[t,:From],
            psd.T[t,:To], psd.MVAbase*sslack_ti[t, 1], psd.MVAbase*sslack_ti[t, 2])
	end    

    return nothing    
end

# method to write the power flow for solve_basecase or solve_contingency

function write_power_flow(OutDir::String, filename::String, psd::SCACOPFdata, p_li::Matrix{Float64}, 
                                p_ti::Matrix{Float64}; cont_idx::Union{Nothing, Int64} = nothing)

    if cont_idx == nothing
        f = open(OutDir * filename, "w")	
        @printf(f, "--base\n")
    elseif cont_idx == 1
        f = open(OutDir * filename, "w")	
        @printf(f, "--contingency\nlabel\n\'%s\'\n", psd.cont_labels[cont_idx])
    else
        f = open(OutDir * filename, "a")	
        @printf(f, "--contingency\nlabel\n\'%s\'\n", psd.cont_labels[cont_idx])
    end
    
	write_power_flow_block(f, psd, p_li, p_ti)
    close(f)

    return nothing
end

# method to write the power flow for solve_SC_ACOPF

function write_power_flow(OutDir::String, filename::String, psd::SCACOPFdata, p_li::Matrix{Float64}, 
                                p_ti::Matrix{Float64}, p_lik::Array{Float64}, p_tik::Array{Float64})

    f = open(OutDir * filename, "w")
    @printf(f, "--base\n")
    write_power_flow_block(f, psd, p_li, p_ti)
    for k = 1:nrow(psd.K)
        @printf(f, "--contingency\nlabel\n\'%s\'\n",
                psd.cont_labels[k])
        write_power_flow_block(f, psd, p_lik[:,:,k], p_tik[:,:,k])
    end
    close(f)

    return nothing
end

# method to write the system slack variable for solve_basecase or solve_contingency

function write_slack(OutDir::String, filename::String, psd::SCACOPFdata, 
                            pslackm_n::Array{Float64}, pslackp_n::Array{Float64}, 
                            qslackm_n::Array{Float64}, qslackp_n::Array{Float64},
                            sslack_li::Matrix{Float64}, sslack_ti::Matrix{Float64};
                            cont_idx::Union{Nothing, Int64} = nothing)

    if cont_idx == nothing
        f = open(OutDir * filename, "w")
        @printf(f, "--base\n")
    elseif cont_idx == 1
        f = open(OutDir * filename, "w")
        @printf(f, "--contingency\n--label\n\'%s\'\n", psd.cont_labels[cont_idx])
    else
        f = open(OutDir * filename, "a")
        @printf(f, "--contingency\n--label\n\'%s\'\n", psd.cont_labels[cont_idx])
    end

    write_slack_block(f, psd, pslackm_n, pslackp_n, qslackm_n, qslackp_n, sslack_li, sslack_ti)
    close(f)

    return nothing
end

# method to write the system slack variable for solve_SC_ACOPF

function write_slack(OutDir::String, filename::String, psd::SCACOPFdata, 
                            pslackm_n::Array{Float64}, pslackp_n::Array{Float64}, 
                            qslackm_n::Array{Float64}, qslackp_n::Array{Float64},
                            sslack_li::Matrix{Float64}, sslack_ti::Matrix{Float64},
                            pslackm_nk::Array{Float64}, pslackp_nk::Array{Float64}, 
                            qslackm_nk::Array{Float64}, qslackp_nk::Array{Float64},
                            sslack_lik::Array{Float64}, sslack_tik::Array{Float64})

    f = open(OutDir * filename, "w")

    @printf(f, "--base\n")
    write_slack_block(f, psd, pslackm_n, pslackp_n, qslackm_n, qslackp_n, sslack_li, sslack_ti)
    for k = 1:nrow(psd.K)
        @printf(f, "--contingency\n--label\n\'%s\'\n",
                psd.cont_labels[k])
        write_slack_block(f, psd, pslackm_nk[:,k], pslackp_nk[:,k], 
                                    qslackm_nk[:,k], qslackp_nk[:,k], 
                                    sslack_lik[:,:,k], sslack_tik[:,:,k])
    end
    close(f)

    return nothing
end

# method to write the components of the power flow constraint

function write_power_flow_cons_block(OutDir::String, filename::String, psd::SCACOPFdata, 
                                sum_pg::Array{Any}, pd::Array{Float64}, gsh::Array{Float64}, 
                                v_n::Array{Float64},  gshvn2::Array{Float64}, sum_p_li::Array{Any},
                                sum_p_ti::Array{Any},  p_relax::Array{Float64}, 
                                sum_qg::Array{Any}, qd::Array{Float64},  Bsh::Array{Float64}, 
                                ssh::Array{Any},  sshvn2::Array{Any}, 
                                sum_q_li::Array{Any}, sum_q_ti::Array{Any}, q_relax::Array{Float64},
                                p_mn::Array{Float64}, p_pn::Array{Float64},
                                q_mn::Array{Float64}, q_pn::Array{Float64};
                                cont_idx::Union{Nothing, Int64} = nothing)

    if cont_idx == nothing
        f = open(OutDir * filename, "w")
        @printf(f, "--base\n")
    elseif cont_idx == 1
        f = open(OutDir * filename, "w")
        @printf(f, "--contingency\n--label\n\'%s\'\n", psd.cont_labels[cont_idx])
    else
        f = open(OutDir * filename, "a")
        @printf(f, "--contingency\n--label\n\'%s\'\n", psd.cont_labels[cont_idx])
    end

    # write active power section
    @printf(f, "--Active power constraint\n")
    @printf(f, "Bus, sum_pg, pd, gsh, v_n, gshvn2, sum_p_li, sum_p_ti, pslackm_n, pslackp_n, p_relax\n")
	for n=1:size(psd.N, 1)
		@printf(f, "%d, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n", 
                psd.N[n,:Bus], sum_pg[n], pd[n], gsh[n], v_n[n], gshvn2[n], sum_p_li[n], sum_p_ti[n], 
                p_mn[n], p_pn[n], p_relax[n])
	end
	
	# write reactive power section
    @printf(f, "--Reactive power constraint\n")
    @printf(f, "Bus, sum_qg, qd, Bsh, Ssh, v_n, BshSshvn2, sum_q_li, sum_q_ti, qslackm_n, qslackp_n, q_relax\n")
	for n=1:size(psd.N, 1)
		@printf(f, "%d, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n", 
                psd.N[n,:Bus], sum_qg[n], qd[n], Bsh[n], ssh[n], v_n[n], sshvn2[n], sum_q_li[n], sum_q_ti[n],
                q_mn[n], q_pn[n], q_relax[n])
	end
    close(f)

    return nothing
end

# function to extract the components of the power flow constraint

function write_power_flow_cons(output_dir::String, filename::String,
                           v_n::AbstractVector, theta_n::AbstractVector,
                           p_li::AbstractArray, q_li::AbstractArray,
                           p_ti::AbstractArray, q_ti::AbstractArray,
                           b_s::AbstractVector,
                           p_g::AbstractVector, q_g::AbstractVector,
                           pslackm_n::AbstractVector, pslackp_n::AbstractVector,
                           qslackm_n::AbstractVector, qslackp_n::AbstractVector,
                           sslack_li::AbstractArray, sslack_ti::AbstractArray,
                           psd::SCACOPFdata, con::GenericContingency=GenericContingency();
                           aux_slack_i::Union{Nothing, AbstractArray}=nothing,
                           cont_idx::Union{Nothing, Int64}=nothing)
    
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
            println("The varible aux_slack_i is not nothing.")
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
            push!(sshvn2, -JuMP.value((-psd.N[n,:Bsh] - sum(b_s[s] for s=psd.SShn[n]))) * JuMP.value(v_n[n]^2))
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

    write_power_flow_cons_block(output_dir, filename, psd, sum_pg, pd, gsh, vn, gshvn2, sum_p_li,
                            sum_p_ti, p_relax, sum_qg, qd, Bsh, ssh, sshvn2, sum_q_li, sum_q_ti, q_relax,
                            p_mn, p_pn, q_mn, q_pn, cont_idx = cont_idx)

    return nothing
end

# method to write the ramp rate information

function write_ramp_rate(OutDir::String, filename::String, psd::SCACOPFdata, r_rate::Array{Float64})

    f = open(OutDir * filename, "w")

    @printf(f, "--ramp rate info\n")
    @printf(f, "Generator, Unit, ramp rate bound\n")
    for g = 1:nrow(psd.G)
		@printf(f, "%d, \'%s\', %.10f\n", psd.G[g,:Bus], psd.G[g,:BusUnitNum], r_rate[g] * psd.MVAbase)
    end
    close(f)

    return nothing
end

# method to write the components of the cost function for solve_SC-ACOPF

function write_cost(OutDir::String, filename::String, psd::SCACOPFdata, bc_gen::Float64, 
                    bc_pen::Float64, con_gen::Array{Float64}, quad_gen::Array{Float64})

	f = open(OutDir * filename, "w")

    @printf(f, "--base cost\n")
    @printf(f, "generation cost, penalty\n")
    @printf(f, "%.10f, %.10f\n", bc_gen, bc_pen)
    @printf(f, "--contingency\n")
    for k = 1:nrow(psd.K)
        @printf(f, "--label\n\'%s\'\n",
                psd.cont_labels[k])
        @printf(f, "contingency penalty, quadratic relaxation \n")
        @printf(f, "%.10f, %.10f\n", con_gen[k], quad_gen[k])
    end
    close(f)

    return nothing
end

# method to write the components of the cost function for solve_basecase

function write_cost(OutDir::String, filename::String, psd::SCACOPFdata, 
                    bc_gen::Float64, bc_pen::Float64)

	f = open(OutDir * filename, "w")

    @printf(f, "--base cost\n")
    @printf(f, "generation cost, penalty\n")
    @printf(f, "%.10f, %.10f\n", bc_gen, bc_pen)
    close(f)
    
    return nothing
end

# method to write the components of the cost function for solve_contingency

function write_cost(OutDir::String, filename::String, psd::SCACOPFdata, 
                    con_gen::Float64, quad_gen::Float64, cont_idx::Int64)

    if cont_idx == 1
        f = open(OutDir * filename, "w")
    else
        f = open(OutDir * filename, "a")
    end

    @printf(f, "--contingency\n")
    @printf(f, "--label\n\'%s\'\n", psd.cont_labels[cont_idx])
    @printf(f, "contingency penalty, quadratic relaxation \n")
    @printf(f, "%.10f, %.10f\n", con_gen, quad_gen)
    close(f)

    return nothing
end
