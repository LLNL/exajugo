# function to write solution block

function write_solution_block(io::IO, psd::SCACOPFdata, sol::SubproblemSolution)
	
	# compute reactive compensaton at each bus
	bcsn = zeros(Float64, nrow(psd.N))
	for ssh = 1:nrow(psd.SSh)
		bcsn[psd.SSh_Nidx[ssh]] += sol.b_s[ssh]
	end
	bcsn *= psd.MVAbase
	
	# write bus section
	@printf(io, "--bus section\n")
    @printf(io, "i, v(p.u.), theta(deg), bcs(MVAR at v = 1 p.u.)\n")
	for n = 1:nrow(psd.N)
		@printf(io, "%d, %.20f, %.20f, %.20f\n", psd.N[n,:Bus],
            sol.v_n[n], 180/pi*sol.theta_n[n], bcsn[n])
	end
	
    # write generator section
	gmap = zeros(Int, nrow(psd.generators))
	for g = 1:nrow(psd.G)
		gmap[psd.G[g,:Generator]] = g
	end
	@printf(io, "--generator section\ni, id, p(MW), q(MW)\n")
	for gi = 1:nrow(psd.generators)
		if gmap[gi] == 0
			@printf(io, "%d, \'%s\', 0, 0\n", psd.generators[gi,:I],
                    psd.generators[gi,:ID])
		else
			g = gmap[gi]
			@printf(io, "%d, \'%s\', %.12f, %.12f\n", psd.G[g,:Bus],
                    psd.G[g,:BusUnitNum], psd.MVAbase*sol.p_g[g],
                    psd.MVAbase*sol.q_g[g])
		end
	end
    
    return nothing
    
end

# method to write base case solution (solution1.txt)

function write_solution(OutDir::String, psd::SCACOPFdata, sol::BasecaseSolution;
                        filename::Union{Nothing, String} = nothing)
    
    # check that solution corresponds to the given data
    if hash(psd) != sol.psd_hash
        error("base case solution does not correspond to power system data.")
    end
    
    # write solution and return
    if filename == nothing
    	f = open(OutDir * "/solution1.txt", "w")
	else
    	f = open(OutDir * filename, "w")
    end
    write_solution_block(f, psd, sol)
	close(f)

    return nothing
    
end

# method to write single contingency solution 

function write_solution(OutDir::String, psd::SCACOPFdata,
                       vsol::ContingencySolution;
                       filename::Union{Nothing, String} = nothing,
                       cont_idx::Int64)
                       
    if cont_idx == 1
        f = open(OutDir * filename, "w")
    else
        f = open(OutDir * filename, "a")
    end

    # write solution and return
    @printf(f, "--contingency\nlabel\n\'%s\'\n",
            psd.cont_labels[cont_idx])
    write_solution_block(f, psd, vsol)
    @printf(f, "--delta section\ndelta(MW)\n%g\n",
            psd.MVAbase*vsol.delta)

    close(f)
    
    return nothing
    
end

# method to write contingency solution (solution2.txt)

function write_solution(OutDir::String, psd::SCACOPFdata,
                       vsol::Vector{ContingencySolution};
                       filename::Union{Nothing, String} = nothing)
    
    # verify that we have a solution for each contingency
    vsol_cont_ids = collect(vsol[i].cont_id for i = 1:length(vsol))
    vsol_cont_ids_unique = sort(vsol_cont_ids)
    unique!(vsol_cont_ids_unique)
    if vsol_cont_ids_unique != psd.K[!,:Contingency]
        missing_cont = setdiff(psd.K[!,:Contingency], vsol_cont_ids_unique)
        @assert length(missing_cont) > 0
        warn_msg = "missing solution data for the following contingencies:"
        for i = 1:length(missing_cont)
            warn_msg *= " "*psd.cont_labels[missing_cont[i]]
        end
        warn_msg *= ". Will proceed only with passed contingencies."
        @warn *
    end
    
    # return here if no contingency solution is passed
    if length(vsol) == 0
        return nothing
    end
    
    # find repeated contingencies if any ...
    perm = sortperm(vsol_cont_ids)
    if length(vsol_cont_ids_unique) == length(vsol)
        unique_vsol_idx = perm
    else
        @warn "found multiple solutions for the same contingency. Will only write first one found."
        unique_vsol_idx = Int[]
        push!(unique_vsol_idx, perm[1])
        for i = 2:length(vsol)
            if vsol_cont_ids[perm[i]] > vsol_cont_ids[perm[i-1]]
                push!(unique_vsol_idx, perm[i])
            else
                @assert vsol_cont_ids[perm[i]] == vsol_cont_ids[perm[i-1]]
            end
        end
    end
    
    # write all passed contingencies
    if filename == nothing
    	f = open(OutDir * "/solution2.txt", "w")
	else
    	f = open(OutDir * filename, "w")
    end
    hash_psd = hash(psd)
    for i = unique_vsol_idx
        if hash_psd != vsol[i].psd_hash
            error("contingency solution does not correspond to SCACOPF data.")
        end
        @printf(f, "--contingency\nlabel\n\'%s\'\n",
                psd.cont_labels[vsol[i].cont_id])
        write_solution_block(f, psd, vsol[i])
        @printf(f, "--delta section\ndelta(MW)\n%g\n",
                psd.MVAbase*vsol[i].delta)
    end
    close(f)
    
    return nothing
    
end

# method to write SCACOPF solution (solution1.txt and solution2.txt)

function write_solution(OutDir::String, psd::SCACOPFdata, sol::SCACOPFsolution;
                        basecase_filename::Union{Nothing, String} = nothing,
                        contingency_filename::Union{Nothing, String} = nothing)
    write_solution(OutDir, psd, sol.basecase, filename = basecase_filename)
    write_solution(OutDir, psd, sol.contingency, filename = contingency_filename)
end
