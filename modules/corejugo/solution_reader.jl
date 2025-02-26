## function to read solution block

function read_solution_block(io::IO)
    
    # check that first line is bus header
    readline(io) == "--bus section" || error("wrong bus section header")
    readline(io);    # ignoring headers
    
    # read bus section
    I_n = Int[]
    v_n = Float64[]
    theta_n = Float64[]
    b_n = Float64[]
    while true
        l = readline(io) 
        l != "--generator section" || break
        fields = split(l, ", ")
        push!(I_n, parse(Int, fields[1]))
        push!(v_n, parse(Float64, fields[2]))
        push!(theta_n, parse(Float64, fields[3]))
        push!(b_n, parse(Float64, fields[4]))
    end
    readline(io);    # ignoring headers
    
    # read generator section
    I_g = Int[]
    ID_g = String[]
    p_g = Float64[]
    q_g = Float64[]
    while true
        l = readline(io)
        (l != "" && l[1:2] != "--") || break
        fields = split(l, ", ")
        push!(I_g, parse(Int, fields[1]))
        push!(ID_g, replace(fields[2], "'" => ""))
        push!(p_g, parse(Float64, fields[3]))
        push!(q_g, parse(Float64, fields[4]))
    end
    
    return I_n, v_n, theta_n, b_n, I_g, ID_g, p_g, q_g
    
end

## function to read base case solution

function read_base_solution(OutDir::String, psd::SCACOPFdata)::BasecaseSolution

    # read base solution file
    io = open(joinpath(OutDir, "solution1.txt"), "r")
    I_n, v_n, theta_n, b_n, I_g, ID_g, p_g, q_g = read_solution_block(io)
    close(io)
    
    # create solution struct to fill in
    sol = BasecaseSolution(psd)
    
    # parse voltage solution
    theta_n .*= pi/180
    Nidx = convert(Vector{Int}, indexin(psd.N[!, :Bus], I_n))
    for n = 1:nrow(psd.N)
        sol.v_n[n] = v_n[Nidx[n]]
        sol.theta_n[n] = theta_n[Nidx[n]]
    end
    
    # parse switched shunt capacitance solution
    b_n ./= psd.MVAbase
    for ssh = 1:nrow(psd.SSh)
        n = psd.SSh_Nidx[ssh]
        if b_n[n] < psd.SSh[ssh, :Blb]
            sol.b_s[ssh] = psd.SSh[ssh, :Blb]
        elseif b_n[n] > psd.SSh[ssh, :Bub]
            sol.b_s[ssh] = psd.SSh[ssh, :Bub]
        else
            sol.b_s[ssh] = b_n[n]
        end
        b_n[n] -= sol.b_s[ssh]
    end
    sum(abs.(b_n)) <= 1E-5 || error("there are ",
                                    round(psd.MVAbase*sum(abs.(b__n)), digits=4),
                                    "MVAR unassigned to shunts.")
    
    # parse generator setpoint solution
    p_g ./= psd.MVAbase
    q_g ./= psd.MVAbase
    Gidx = convert(Vector{Int},
                   indexin(string.(psd.G[!, :Bus], ":", psd.G[!, :BusUnitNum]),
                           string.(I_g, ":", ID_g)))
    for g = 1:nrow(psd.G)
        sol.p_g[g] = p_g[Gidx[g]]
        sol.q_g[g] = q_g[Gidx[g]]
    end
    
    # return solution struct
    return sol
    
end
