#__precompile__()

module CompetitionInterface

## elements to be exported

export MyJulia1_main, MyJulia2_main

## load external modules
using Printf
using Random, Ipopt, JuMP, DataFrames
import DataStructures: OrderedDict

using MPI#, Distributed

## load internal modules
Base.include(@__MODULE__,"InstanceReader.jl")
Base.include(@__MODULE__,"SmoothApproximations.jl")
Base.include(@__MODULE__,"SCACOPF.jl")
Base.include(@__MODULE__,"SolutionEvaluator.jl")
Base.include(@__MODULE__,"SolutionWriter.jl")
Base.include(@__MODULE__,"GoUtils.jl")
Base.include(@__MODULE__,"SolutionReader.jl")
using .InstanceReader, .SmoothApproximations, .SCACOPF, .SolutionEvaluator, .SolutionWriter, .GoUtils, .SolutionReader

## constants
import .SmoothApproximations: DEFAULTMU

##
function selectContingencies(NumCont::Int, G::DataFrame, K::DataFrame)
    conidx = []
    if NumCont>4 || size(K,1)<=4
        l = findfirst(K[:ConType] .== :Line)
        if l!=nothing; append!(conidx, l); end
        t = findfirst(K[:ConType] .== :Transformer)
        if t!=nothing; append!(conidx, t); end
    end

    ngencont = NumCont-size(conidx,1)

    if ngencont>0
        # indexes of generator conting
        genconidx = findall(K[:ConType] .== :Generator)
        if !isempty(genconidx)
            # IDs of out generators
            genoutid = K[:IDout][genconidx]
            
            # indexes in G of out generators 
            genoutidx = findall((in)(genoutid), G[:Generator])
            
            # indexes (in G[out]) of largest out generators
            genlargeidx = partialsortperm(G[:Pub][genoutidx], 1:ngencont, rev=true)
            
            # indexes in G of largest out generators 
            genoutlargeidx = genoutidx[genlargeidx]
            
            contgenlargeidx=findall((in)(G[:Generator][genoutlargeidx]), genoutid)
            
            @assert size(contgenlargeidx,1) == ngencont
            
            append!(conidx, genconidx[contgenlargeidx])
        end
    end
    return conidx
end

function writeSolution1(myRank, OutDir, SCACOPFSol, MVAbase, N, SSh, G, generators)
    # write base case solution ("solution1.txt", using GO competition format)
    clock = @elapsed begin
	mkpath(OutDir)
        writesolution(OutDir, 
                      SCACOPFSol.base_state.v_n, SCACOPFSol.base_state.theta_n,
                      SCACOPFSol.base_state.b_s, SCACOPFSol.base_state.p_g,
                      SCACOPFSol.base_state.q_g, MVAbase, N, SSh, G, generators)
        
    end
    println("Base case solution written in ", round(clock, digits=2), " secs on rank ", myRank)
    return nothing
end

# priority of writing the solutions is given by the rank in the communicator
#
# the problems assigned to each rank are such that 
# - lower ranks are likely to finish sooner than higher ranks
# - rank 0 has a tight max num of iterations and solves the smallest 
# problem, which is (needs to be) chosen to finish with certainty
# Parameters:
# - boptimal: indicates whether local 'myRank' rank has a valid solution
function comm_getMaxPriority(myRank, boptimal)
    commSz = MPI.Comm_size(MPI.COMM_WORLD)
    nWorkers = commSz-1
    if myRank==0
        #master
        workersOnline = 1:nWorkers
        
        #we assume (precondition) that rank 0 has already written solution1
        maxPriority = zeros(Int64, 1)
        workersRecvBuf = zeros(Int64, nWorkers)
        requests = Array{MPI.Request}(undef, nWorkers)
        for r in workersOnline
            requests[r] = MPI.Irecv!(view(workersRecvBuf, r), 1, r, 117+r, MPI.COMM_WORLD)
        end
        
        while length(workersOnline)>0
            sleep(0.05)
            workersDone = []
            for r in workersOnline
                (done,reqStat) = MPI.Test!(requests[r])
                if done
                    #send maxPriority rank that already written solution1.txt; 
                    #workers with lower priority (or that did not solve to optimality) will not write
                    MPI.Send(maxPriority, 1, r, 10117+r, MPI.COMM_WORLD)
                    
                    if workersRecvBuf[r]==1
                        maxPriority[1] = max(maxPriority[1], r)
                    else
                        @printf("warning: recv on master not 1 - worker %d did not solve to optimality", r)
                    end
                    sleep(0.01)
                    push!(workersDone, r)
                end
            end
            
        workersOnline = setdiff(workersOnline, workersDone)
        end
        
        @printf("Master rank %d finished with maxPriorityCompleted=%d\n", myRank, maxPriority[1])
        return maxPriority[1]
    else
        #workers
        completed = ones(Int64,1)
        if !boptimal; completed[1]=0; end

        request = MPI.Isend(completed, 1, 0, 117+myRank, MPI.COMM_WORLD)
        reqStat = MPI.Wait!(request)

        #@show stat
        maxPriorityCompleted=ones(Int64,1)

        recvStat = MPI.Recv!(maxPriorityCompleted, 1, 0, 10117+myRank, MPI.COMM_WORLD)
        @printf("Work rank %d finished with maxPriorityCompleted=%d\n", myRank, maxPriorityCompleted[1])
        return maxPriorityCompleted[1]
    end
end

## MyJulia1 function as defined in: https://gocompetition.energy.gov/languages
function MyJulia1_main(InFile1::String, InFile2::String, InFile3::String, InFile4::String,
	          TimeLimitInSeconds::Real, ScoringMethod::Int, NetworkModel::String, 
                  tmstart::UInt64)
    println("MyJulia1 called with:\n",
	    "\tInFile1 .......... : ", InFile1, "\n",
	    "\tInFile2 .......... : ", InFile2, "\n",
	    "\tInFile3 .......... : ", InFile3, "\n",
	    "\tInFile4 .......... : ", InFile4, "\n",
	    "\tTimeLimitInSeconds : ", TimeLimitInSeconds, "\n",
	    "\tScoringMethod .... : ", ScoringMethod, "\n",
	    "\tNetworkModel ..... : ", NetworkModel)
    o = GoOptions()
    o.PostSolveContingencies = false
    o.tmstart = tmstart
    io = stdout # set up log file
    OutDir = "./"

    flush(io)

    ENV["OMP_NUM_THREADS"]=1
    ENV["MKL_NUM_THREADS"]=1

    MPI.Init()
    myRank = MPI.Comm_rank(MPI.COMM_WORLD)
    commSz = MPI.Comm_size(MPI.COMM_WORLD)
    if myRank==0
        omp = haskey(ENV,"OMP_NUM_THREADS") ? ENV["OMP_NUM_THREADS"] : "NA"
        mkl = haskey(ENV,"MKL_NUM_THREADS") ? ENV["MKL_NUM_THREADS"] : "NA"
        @printf("A run involving %d ranks: OMP_NUM_THREADS=%s MKL_NUM_THREADS=%s\n", commSz, omp, mkl)
    end

    # read instance
    println(io, "Reading and parsing instance data ..."); flush(io)
    clock = @elapsed begin
	MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches,
	switchedshunts, generatordsp, activedsptables, costcurves, governorresponse,
	contingencies = readinstance(InFile3, InFile4, InFile2, InFile1)
	N, L, T, SSh, G, K, P = 
            GOfmt2params(MVAbase, buses, loads, fixedbusshunts, generators,
			 ntbranches, tbranches, switchedshunts, generatordsp, activedsptables, costcurves,
			 governorresponse, contingencies)
    end
    if myRank==0
        println(io, "Data prepared for problem formulation in ", round(clock, digits=2), " secs.")
        flush(io)
    end

    if myRank==0
        @printf("Generators %d   Lines %d   Buses %d  Contingencies %d\n", size(G,1), size(L,1), size(N,1), size(K,1))
    end

    # compute index and sets for performant formulation
    L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn, K_outidx = indexsets(N, L, T, SSh, G, K)

    o.CouplingMode = :complementarityapprox #:fischerburm
    o.SmoothingParamCoupling=1e-2

    rankActive=false

    ## common to all ranks and problem sizes
    o.ipopt.linear_solver = "ma57"
    o.ipopt.mu_init = 1
    o.ipopt.tol           = 1e-7
    o.ipopt.mu_target     = 5e-9
    o.ipopt.print_level   = 5
    o.ipopt.print_frequency_iter=10

    ##common to all problem sizes
    if myRank==0
        o.ipopt.max_iter = 400
        o.CouplingMode = :ignore
        rankActive = true
    elseif myRank==1
        o.CouplingMode = :complementarityapprox
        o.SmoothingParamCoupling = 1e-2
        rankActive = true
    elseif myRank==2
        o.CouplingMode = :complementarityapprox
        o.SmoothingParamCoupling = 1e-2
        rankActive = true
    elseif myRank==3
        o.CouplingMode = :complementarityapprox
        o.SmoothingParamCoupling = 1e-3
        rankActive = true
    else #5 or more ranks
        #not active
    end

    ## sizes less than 2K
    NumCont = (ScoringMethod==1 || ScoringMethod==3) ? 10 : 16;
    if myRank==1
        #same NumCont as rank 0
        # with smoothing
    elseif myRank==2
        NumCont += (ScoringMethod==1 || ScoringMethod==3) ? 2 : 4
    elseif myRank==3
        NumCont += (ScoringMethod==1 || ScoringMethod==3) ? 2 : 4
        #smoothing is tighter at 1e-3
    else #5 or more ranks
        #not active
    end
        
    ## sizes between 2K and 5K
    if size(L,1)>1999 || size(N,1)>1999
        NumCont = (ScoringMethod==1 || ScoringMethod==3) ? 2 : 8
        if myRank==1
            # same NumCont as rank 0
            # but with smoothing
        elseif myRank==2
            NumCont += (ScoringMethod==1 || ScoringMethod==3) ? 2 : 4
        elseif myRank==3
            NumCont += (ScoringMethod==1 || ScoringMethod==3) ? 2 : 4
            #smoothing is tighter at 1e-3
        else #5 or more ranks
            #not active
        end
    end

    if size(L,1)>4999 || size(N,1)>4999
        NumCont = (ScoringMethod==1 || ScoringMethod==3) ? 1 : 5;
       if myRank==1
            # same NumCont as rank 0
            # but with smoothing
        elseif myRank==2
            NumCont += (ScoringMethod==1 || ScoringMethod==3) ? 2 : 3
        elseif myRank==3
            NumCont += (ScoringMethod==1 || ScoringMethod==3) ? 2 : 3
            #smoothing is tighter at 1e-3
        else #5 or more ranks
            #not active
        end
    end

    if size(L,1)>9999 || size(N,1)>9999
        NumCont = (ScoringMethod==1 || ScoringMethod==3) ? 1 : 3;
       if myRank==1
            # same NumCont as rank 0
            # but with smoothing
        elseif myRank==2
            NumCont += (ScoringMethod==1 || ScoringMethod==3) ? 1 : 2
        elseif myRank==3
            NumCont += (ScoringMethod==1 || ScoringMethod==3) ? 1 : 2
            #smoothing is tighter at 1e-3
        else #5 or more ranks
            #not active
        end
    end
    
    if size(L,1)>14999 || size(N,1)>14999
        NumCont = (ScoringMethod==1 || ScoringMethod==3) ? 0 : 1;
        if myRank==1
            NumCont = (ScoringMethod==1 || ScoringMethod==3) ? 1 : 2;
            # but with smoothing
       elseif myRank==2
            NumCont += (ScoringMethod==1 || ScoringMethod==3) ? 2 : 3
        elseif myRank==3
            NumCont += (ScoringMethod==1 || ScoringMethod==3) ? 2 : 3
            #smoothing is tighter at 1e-3
        else #5 or more ranks
            #not active
       end
    end

    boptim = false

    if rankActive
        # select up to NumCont contingencies in conidx 
        NumCont = min(NumCont, size(K,1))    
        conidx = selectContingencies(NumCont, G, K)
        
        SCACOPFSol = nothing
        
        sleep(0.01*myRank)
        print(io, "On rank ", myRank, " contingencies considered: ", K[conidx,:], "\n"); flush(io)
        
        clock1 = @elapsed begin       
            
            tm_since_start = (time_ns() - o.tmstart)/1e9

            o.ipopt.max_cpu_time  = TimeLimitInSeconds-tm_since_start-4

            #@show o.ipopt.max_cpu_time
            #@show tm_since_start

            SCACOPFSolver = with_ipopt_optimizer(o)
        
            try
	        SCACOPFSol1 = SolveSCACOPF(o, N, L, T, SSh, G, K[conidx,:], P, 
                                           SCACOPFSolver,
			                   IndexSets=nothing)
                SCACOPFSol = SCACOPFSol1
                boptim = true
            catch ex
                println(ex)
                println("exception in pass 1 on rank ", myRank); flush(io)
            end
        end
        println(io, "Pass 1 problem solved in ", round(clock1, digits=2), " secs on rank ", myRank); flush(io)
        tm_since_start = (time_ns() - o.tmstart)/1e9
        @show tm_since_start
        flush(io)
        if myRank==0
            if boptim
                #master always right the solution before entering the communication as a backup
                writeSolution1(myRank, OutDir, SCACOPFSol, MVAbase, N, SSh, G, generators)
            else
                println(io, "rank 0: could not write solution - solve failed"); flush(io)
            end
        end
    end # of rankActive==true

    ## communication among all ranks in COMM_WORLD (!including non-active ones)
    writtenSolPriority = comm_getMaxPriority(myRank, boptim)
    if boptim && myRank>writtenSolPriority
        writeSolution1(myRank, OutDir, SCACOPFSol, MVAbase, N, SSh, G, generators)
    end
    

    MPI.Barrier(MPI.COMM_WORLD)
    sleep(0.001*myRank)
    println(io, "MyJulia1 took ", (time_ns()-o.tmstart)/1e9, " sec on rank ", myRank); flush(io)
    MPI.Finalize()
end

## MyJulia2 function as defined in: https://gocompetition.energy.gov/languages


function MyJulia2_main(InFile1::String, InFile2::String, InFile3::String, InFile4::String,
	          TimeLimitInSeconds::Real, ScoringMethod::Int, NetworkModel::String, 
                  tmstart::UInt64; sol1Dir::String=".")
    println("MyJulia2 called with:\n",
	    "\tInFile1 .......... : ", InFile1, "\n",
	    "\tInFile2 .......... : ", InFile2, "\n",
	    "\tInFile3 .......... : ", InFile3, "\n",
	    "\tInFile4 .......... : ", InFile4, "\n",
	    "\tTimeLimitInSeconds : ", TimeLimitInSeconds, "\n",
	    "\tScoringMethod .... : ", ScoringMethod, "\n",
	    "\tNetworkModel ..... : ", NetworkModel)

    io = stdout # set up log file
    OutDir = "./"; flush(io)

    ENV["OMP_NUM_THREADS"]=1
    ENV["MKL_NUM_THREADS"]=1

    MPI.Init()
    myRank = MPI.Comm_rank(MPI.COMM_WORLD)
    commSz = MPI.Comm_size(MPI.COMM_WORLD)

    if(commSz==1)
        println("Need at least two ranks to run MJ2. Got 1.")
        return
    end

    if myRank==0
        omp = haskey(ENV,"OMP_NUM_THREADS") ? ENV["OMP_NUM_THREADS"] : "NA"
        mkl = haskey(ENV,"MKL_NUM_THREADS") ? ENV["MKL_NUM_THREADS"] : "NA"
        @printf("A run involving %d ranks: OMP_NUM_THREADS=%s MKL_NUM_THREADS=%s\n", commSz, omp, mkl); flush(io)
    end    
    
    # read instance
    println(io, "Reading and parsing instance data ..."); flush(io)
    clock = @elapsed begin
	MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches,
	switchedshunts, generatordsp, activedsptables, costcurves, governorresponse,
	contingencies = readinstance(InFile3, InFile4, InFile2, InFile1)

	N, L, T, SSh, G, K, P = 
            GOfmt2params(MVAbase, buses, loads, fixedbusshunts, generators,
			 ntbranches, tbranches, switchedshunts, generatordsp, activedsptables, costcurves,
			 governorresponse, contingencies)
    end
    if myRank==0; @printf("Generators %d   Lines %d   Buses %d  Contingencies %d\n", size(G,1), size(L,1), size(N,1), size(K,1)) end
    println(io, "Data prepared for problem formulation in ", round(clock, digits=2), " secs on rank ", myRank); flush(io)



    # compute index and sets for performant formulation
    L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn, K_outidx = indexsets(N, L, T, SSh, G, K)
    
    # read base case solution
    println(io, "Reading base case solution from 'solution1.txt' ...")
    clock = @elapsed begin
	v_n, theta_n, b_s, p_g, q_g = ReadBaseSolution(sol1Dir, MVAbase, N, SSh, G)
	_, prodcost, penaltybase = BaseCaseValue(v_n, theta_n, b_s, p_g, q_g,
					      N, L, T, SSh, G, P)
	SCACOPFSol = SCACOPF.SCACOPFState(SCACOPF.BaseState(v_n, theta_n, b_s,
							    p_g, q_g, prodcost,
							    penaltybase))
    end
    println(io, "Base case solution read in ", round(clock, digits=2), " secs on rank ", myRank, " globtime ", (time_ns()-tmstart)/1e9); flush(io)
    
    # initial allocation of contingencies for workers/slaves
    #num contingencies in each batch sent to ranks (for rank 0 there are some exceptions)


    nK = size(K, 1)
    contingencies_left = collect(1:nK) #only on rank 0
    contingPerRank = div(nK,commSz-1);
    contingPerRank = max(contingPerRank,1);

    sz_con_batch=min(1, contingPerRank)

    q_con_local = Vector{Int64}()
    con_slaves= Dict{Int64,Vector{Int64}}() #this is for rank 0

    if myRank>0
        for c=1:sz_con_batch
            i = 1  #everybody solve conting 1 during warmup
            #i = (myRank-1)*contingPerRank + c
            if i<=nK
                push!(q_con_local, i)
            end
        end
    else
        #nothing master does not solve
    end

    #mapping inside the master (myRank==0)
    if myRank==0
        con_slaves[0]=Vector{Int64}() #empty for master
        for r=1:commSz-1
            con_slave = (r-1)*contingPerRank .+ collect(1:sz_con_batch)
            con_slave = con_slave[con_slave .<= nK]
            con_slaves[r] = con_slave
            for c in con_slave
                deleteat!(contingencies_left, findall(isequal(c), contingencies_left))
            end
        end
    end

    if myRank==0
        @show con_slaves
    end
@show (myRank, q_con_local)

    StartingPoint = nothing
    ContSol = nothing
    #all workers solve the first assigned problem to get a StartingPoint and sizes
    pass1time = @elapsed begin
        k = nothing
        if size(q_con_local,1)>0
            k = q_con_local[1]
        elseif myRank==0 
            #force warming with 1 contingency on rank 0 -> needed to get the sizes
            k = 1
        end
        if k!=nothing
            # define options and other auxiliary vars
	    o = SCACOPF.GoOptions()
	    o.tmstart = tmstart
	    o.CouplingMode = :complementarityapprox
	    o.SmoothingParamCoupling=1e-2
	    o.ipopt.mu_init = 1
	    o.ipopt.linear_solver = "ma57"
	    o.ipopt.tol           = 1e-8
	    o.ipopt.mu_target     = 1e-9
	    o.ipopt.print_level = 0
            o.ipopt.max_cpu_time = 5000.
            o.ipopt.max_iter = 1500
	    SolverCon = SCACOPF.with_ipopt_optimizer(o)
            
            println("[warmup] contingency ", k, " out of ", size(K, 1), " on rank ", myRank); flush(io)

	    contime = @elapsed begin
		ContSol, StartingPoint = 
                    SCACOPF.SolveContingencyTwoStep(o, N, L,
				                    T, SSh, G, (K[:ConType][k], K[:IDout][k], K[:Contingency][k]), P,
				                    SCACOPFSol, SolverCon, StartingPoint=StartingPoint,
				                    IndexSets=(L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn,
				                               Lin, Tidxn, Tin, SShn, Gn, K_outidx[k]))
		@assert ContSol.cont_id == K[:Contingency][k]
		end
		println("[warmup] contingency ", k, " out of ",  size(K, 1), " solved in ",
			round(contime, digits=3), " secs on ", myRank); flush(io)
        else
            println("Rank ", myRank, " does not have contingencies assigned"); flush(io)
        end

    end # of pass1time = @elapsed begin
    println("Warming up on rank ", myRank, " took ", round(pass1time, digits=3), " secs ", " globtime ", (time_ns()-tmstart)/1e9); flush(io)

    #reset local contingency list
    q_con_local = Vector{Int64}()
    if myRank>0
        for c=1:sz_con_batch
            i = (myRank-1)*contingPerRank + c
            if i<=nK
                push!(q_con_local, i)
            end
        end
    else
        #nothing master does not solve
    end

    nDouComm = ContSol==nothing ? 0 : SCACOPF.getTotalSize(ContSol)

    # master creates Irecv requests
    req_contsol_recv = Dict{Int64, Vector{MPI.Request}}()
    contsol_buffer = []
    if myRank==0
        contsol_buffer = Dict{Int64, Vector{Float64}}();
        for rank=0:commSz-1
            req_contsol_recv[rank] = Vector{MPI.Request}()
            
            for i=1:size(con_slaves[rank],1)

                c = con_slaves[rank][i]
                contsol_buffer[c] = zeros(Float64, nDouComm)
                
                tag2 = 2 ## =c
                request = MPI.Irecv!(contsol_buffer[c], nDouComm, rank, tag2, MPI.COMM_WORLD)
                @printf("[rank 0] contsol recv created for conting %d on rank %d\n", c, rank)
                push!(req_contsol_recv[rank], request)
            end
        end
    end

   # myRank>0 solves the contingencies and (I)send them as soon they're finished
   # everybody requests more contingencies when q_con_local has 2 or less contingencies
   contsol_buf_send = Dict{Int64, Vector{Float64}}()

   req_contsol_send = Vector{MPI.Request}()
   #aaa

   req_contidxs_send = Vector{MPI.Request}()
   req_contidxs_recv = Vector{MPI.Request}()
   contidxs_buf_send = zeros(Int64, 1)
   contidxs_buf_recv = zeros(Int64, sz_con_batch)

   # only rank 0 uses these
   req_contidxs_recv0 = Dict{Int64, Vector{MPI.Request}}()
   contidxs_buf_recv0 = Dict{Int64, Vector{Int64}}()
   req_contidxs_send0 = Dict{Int64, Vector{MPI.Request}}()
   contidxs_buf_send0 = Dict{Int64, Vector{Int64}}()

   if myRank==0
     for r=0:commSz-1
         req_contidxs_recv0[r] = Vector{MPI.Request}()
         req_contidxs_send0[r] = Vector{MPI.Request}()
         contidxs_buf_recv0[r] = zeros(Int64, 1)
         contidxs_buf_send0[r] = zeros(Int64, sz_con_batch)
     end
   end

   comm_contidxs_done = zeros(Int64, commSz);

   cont_req_num = 1
   cont_req_num0= ones(Int64,commSz)

   ask_for_conting=true
   allcontime = @elapsed begin
       while true

           if size(q_con_local,1)>0
@show myRank,q_con_local
               k = q_con_local[1]
               deleteat!(q_con_local,1)
               
               # define options and other auxiliary vars
               o = SCACOPF.GoOptions()
               o.tmstart = tmstart
               o.CouplingMode = :complementarityapprox
               o.SmoothingParamCoupling=1e-4
               o.ipopt.mu_init = 1
               o.ipopt.linear_solver = "ma57"
               o.ipopt.tol           = 1e-8
               o.ipopt.mu_target     = 1e-9
               o.ipopt.print_level = 0
               o.ipopt.max_cpu_time = 3600.
               o.ipopt.max_iter = 1000
               SolverCon = SCACOPF.with_ipopt_optimizer(o)
               
               println("Solving contingency ", k, " out of ", size(K, 1), " on rank ", myRank, "  globtime ", (time_ns()-tmstart)/1e9); flush(io)

               contime = @elapsed begin
	           ContSol, StartingPoint = 
                       SCACOPF.SolveContingencyTwoStep(o, N, L,
				                       T, SSh, G, (K[:ConType][k], K[:IDout][k], K[:Contingency][k]), P,
				                       SCACOPFSol, SolverCon, StartingPoint=StartingPoint,
				                       IndexSets=(L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn,
				                                  Lin, Tidxn, Tin, SShn, Gn, K_outidx[k]))
	           @assert ContSol.cont_id == K[:Contingency][k]
               end
               println("\nContingency ", k, " out of ",  size(K, 1), " solved in ",
	               round(contime, digits=3), " secs on rank ", myRank, " globtime ", (time_ns()-tmstart)/1e9); flush(io)
               
               # send ContSol
               contsol_buf_send[k] = SCACOPF.toArray(ContSol);
               @assert size(contsol_buf_send[k], 1) == nDouComm
               
               tag2 = 2 ## =k
               request  = MPI.Isend(contsol_buf_send[k], size(contsol_buf_send[k],1), 0, tag2, MPI.COMM_WORLD)
               push!(req_contsol_send, request)
           end #size(q_con_local,1)>0

           # master part -> after each solve, master checks on communication
           if myRank==0
               for rank=0:commSz-1
                   for i=1:size(con_slaves[rank],1)
                   
                       c = con_slaves[rank][i]

                       #check if already completed
                       if c<0; continue; end
                   
                       (recvdone,reqStat) = MPI.Test!(req_contsol_recv[rank][i])
                       if recvdone
                           @printf("[on rank 0] contsol recv from rank=%d conting=%d done=%d\n", rank, c, recvdone);
                           con_slaves[rank][i]=-1           
                       end
                   end
               end
               flush(io)
           end #if myRank==0

           # workers part - if q_con_local is small, request additional contingencies -> msg id nK+1 +myRank*nK+ ...
           if myRank>0
               if size(q_con_local,1)<=2 && ask_for_conting 
                   #there may be data in the pipe - only if the send and recv completed
                   if isempty(req_contidxs_recv) && isempty(req_contidxs_send)
                       contidxs_buf_send[1] = cont_req_num
                       
                       tag3 = 3333 # = nK+1+myRank*nK+cont_req_num
                       request = MPI.Isend(contidxs_buf_send, 1, 0, tag3, MPI.COMM_WORLD)
                       @printf("on rank %d: send request %d for contingencies globtime %g\n", myRank, cont_req_num, (time_ns()-tmstart)/1e9); flush(io)
                       push!(req_contidxs_send, request)
                   end
               end 
           end

           # workers part -> check on requests for contidxs -> msg id nK+1 +comSz*nK + myRank*nK+ ...
           if myRank>0
               if ask_for_conting && !isempty(req_contidxs_send)
                   (senddone,reqStat) = MPI.Test!(req_contidxs_send[end])
                   if senddone
                       @printf("on rank %d: send request %d for contingencies completed globtime %g\n", myRank, cont_req_num, (time_ns()-tmstart)/1e9)
                   
                       # post the recv for the actual indexes
                       msg_id = cont_req_num 

                       tag4 = 4444 # = nK+1+ commSz*nK + myRank*nK+msg_id
                       request = MPI.Irecv!(contidxs_buf_recv, sz_con_batch, 0, tag4, MPI.COMM_WORLD)
                       push!(req_contidxs_recv, request)
                       @printf("on rank %d: irecv request %d for indexes globtime %g\n", myRank, msg_id, (time_ns()-tmstart)/1e9)
                       @assert size(req_contidxs_recv,1) == size(req_contidxs_send,1)
                       @assert size(req_contidxs_recv,1) == 1
                       
                       deleteat!(req_contidxs_send, size(req_contidxs_send,1))
                   
                       flush(io)
                   end
               end #if ask_for_conting && !isempty(req_contidxs_send)
           end 

           if !isempty(req_contidxs_recv)
               (recvdone, reqstat) = MPI.Test!(req_contidxs_recv[end])
               if recvdone
                   @printf("on rank %d: irecv request %d for indexes completed globtime %g\n", myRank, cont_req_num, (time_ns()-tmstart)/1e9)

                   for cc in contidxs_buf_recv
                       if cc>0
                           push!(q_con_local, cc)
                       end
                       if cc==-1
                           ask_for_conting=false
                       end
                   end
                   deleteat!(req_contidxs_recv, size(req_contidxs_recv,1))
                   cont_req_num += 1
               end
           end

           
           # master part - 
           # i. check for requests for additional contingencies; 
           # ii. (I)send these contingencies
           # iii. create contsol recv requests for the contingencies that were just sent 
           if myRank==0
               all_comm_contidxs_done=true
               for r=commSz-1:-1:1
                   if 1==comm_contidxs_done[r+1]; continue; end

                   all_comm_contidxs_done=false

                   if isempty(req_contidxs_recv0[r]) && isempty(req_contidxs_send0[r])
                       # initiate Irecv 
                       contidxs_buf_recv0[r][1]=0
                       msg_id = cont_req_num0[1+r]
                       @printf("[rank 0] recv request %d for contingencies for rank %d\n", msg_id, r)
                       tag3  = 3333 # = nK+1+r*nK+msg_id
                       request = MPI.Irecv!(contidxs_buf_recv0[r], 1, r, tag3, MPI.COMM_WORLD)
                       push!(req_contidxs_recv0[r], request)
                   end

                   #test  contidxs_buf_recv[r]    
                   if !isempty(req_contidxs_recv0[r])
                       (recvdone,reqStat) =  MPI.Test!(req_contidxs_recv0[r][end])
                       if recvdone
                           msg_id = cont_req_num0[1+r]
                           @printf("[rank 0] recv request %d for contingencies for rank %d completed\n", msg_id, r)
                           
                           contidxs_buf_send0[r] = -ones(Int64, sz_con_batch)
                           ncont = min(sz_con_batch, size(contingencies_left,1))
                           
                           selection = contingencies_left[ contingencies_left .>= (r-1)*contingPerRank+1 ]
                           if size(selection,1)<ncont
                               #look from the start
                               selection = contingencies_left[ contingencies_left .>= 0 ]
                           end
                           
                           contidxs_buf_send0[r][1:ncont] = selection[1:ncont]
                           for c in contidxs_buf_send0[r][1:ncont]
                               deleteat!(contingencies_left, findall(isequal(c), contingencies_left))
                           end
                           
                           # (ii.)
                           msg_id = cont_req_num0[1+r]

                           tag4 =  4444 # = nK+1+commSz*nK+r*nK+msg_id
                           request = MPI.Isend(contidxs_buf_send0[r], sz_con_batch, r, tag4, MPI.COMM_WORLD)
                           @printf("[rank 0] send request %d for indexes on rank %d\n", msg_id, r)
                           push!(req_contidxs_send0[r], request)
                           @assert size(req_contidxs_recv0[r],1) == size(req_contidxs_send0[r],1)
                           
                           #delete recv request
                           deleteat!(req_contidxs_recv0[r], size(req_contidxs_recv0[r],1))
                           
                           #
                           #create contsol recvs for the contingencies that were just sent 
                           #            
                           for c in contidxs_buf_send0[r] 
                               if c<=0; continue; end
                               push!(con_slaves[r],c)
                               contsol_buffer[c] = zeros(Float64, nDouComm)
                               
                               tag2 = 2 ## = c
                               request = MPI.Irecv!(contsol_buffer[c], nDouComm, r, tag2, MPI.COMM_WORLD)
                               @printf("[rank 0] contsol recv created for conting %d on rank %d\n", c, r)
                               push!(req_contsol_recv[r], request)
                           end
#@show con_slaves
                           if ncont<sz_con_batch
                               #apparently we are out of contingencies
                               @assert 0==size(contingencies_left,1)
                               comm_contidxs_done[r+1]=1
                           end
                           
                       end #if recvdone
                   end #if !isempty(req_contidxs_recv0[r])

                   # test req_contidxs_send0[r]
                   if !isempty(req_contidxs_send0[r])
                       (senddone, reqStat) = MPI.Test!(req_contidxs_send0[r][end])
                       if senddone
                           msg_id = cont_req_num0[1+r]
                           @printf("[rank 0] send request %d for indexes on rank %d completed\n", msg_id, r)
                           cont_req_num0[1+r] += 1
                           deleteat!(req_contidxs_send0[r], size(req_contidxs_send0[r],1))
                       end
                   end #if !isempty(contidxs_send0[r])
               end

               # myRank==0
               if all_comm_contidxs_done && size(q_con_local,1)==0 && isempty(req_contidxs_recv) && isempty(req_contidxs_send) 
                   sleep(0.001)
                   flush(io)
                   break
               end
               flush(io)
           else # myRank>0
               if !ask_for_conting && size(q_con_local,1)==0 && isempty(req_contidxs_recv) && isempty(req_contidxs_send) 
                   #sleep(0.5)
                   flush(io)
                   break
                end
               flush(io)
           end #if myRank==0


           if size(q_con_local,1)==0
               sleep(0.001)
           end
       end # of while
  end # of allcontime = @elapsed begin
  println(io, "All contingencies solved in ", round(allcontime, digits=2), " secs on rank ", myRank); flush(io)
 
  # map solution vector into ordered dictionary (this only creates pointers)
  if myRank == 0
      ContSolutions = OrderedDict{Int, SCACOPF.ContingencyState}()
      
      # wait for all to complete
      done = false
      while !done
          done = true
          for rank=0:commSz-1
              for i=1:size(con_slaves[rank],1)
                  
                  c = con_slaves[rank][i]
                  
                  #check if already completed
                  if c<0; continue; end
                  
                  (recvdone,reqStat) = MPI.Test!(req_contsol_recv[rank][i]) 
                  if recvdone
                      @printf("[on rank 0] contsol recv from rank=%d conting=%d done=%d\n", rank, c, recvdone);
                      con_slaves[rank][i]=-1                           
                  else
                      @printf("[on rank 0] contsol recv from rank=%d conting=%d NOT done\n", rank, c);
                      done=false;
                  end
              end
          end
          if !done; sleep(0.1); end
      end
      println("Master: all recv's are completed");
      for k = 1:size(K, 1)
          # ContSol is used as a template
          ContSolutions[K[:Contingency][k]] = SCACOPF.fromArray(contsol_buffer[k], ContSol)
          @assert ContSolutions[K[:Contingency][k]].cont_id == K[:Contingency][k]
          contsol_buffer[k] = []
      end
      
      
      # write contingency solution ("solution2.txt", using GO competition format)
      println(io, "Writting contingencies solutions ...")
      writetime = @elapsed begin
          #! for now just copy sol to the arrays so that the IO code is not touched
          #! we will revisit this later
          
          @assert(length(ContSolutions)==size(K,1))
          
	  v_nk     = Array{Float64, 2}(undef, size(N, 1), size(K, 1))
	  theta_nk = Array{Float64, 2}(undef, size(N, 1), size(K, 1))
	  b_sk     = Array{Float64, 2}(undef, size(SSh, 1), size(K, 1))
	  p_gk     = Array{Float64, 2}(undef, size(G, 1), size(K, 1))
	  q_gk     = Array{Float64, 2}(undef, size(G, 1), size(K, 1))
	  delta_k  = Vector{Float64}(undef, size(K, 1))
          
          for k=1:size(K,1)
              ContSol = ContSolutions[K[:Contingency][k]]
              @assert K[:Contingency][k]==ContSol.cont_id
              
              v_nk[:,k]     = ContSol.v_n
              theta_nk[:,k] = ContSol.theta_n
              b_sk[:,k]     = ContSol.b_s
              p_gk[:,k]     = ContSol.p_g
              q_gk[:,k]     = ContSol.q_g
              delta_k[k]    = ContSol.delta
          end
          
	  writesolution(OutDir,
                        SCACOPFSol.base_state.v_n, 
                        SCACOPFSol.base_state.theta_n,
                        SCACOPFSol.base_state.b_s,
                        SCACOPFSol.base_state.p_g,
                        SCACOPFSol.base_state.q_g,  
                         MVAbase, N, SSh, G, generators, true, 
                        v_nk, theta_nk, b_sk,
		        p_gk, q_gk, delta_k, K, contingencies)
      end
      println(io, "Contingencies solutions written in ", round(writetime, digits=2), " secs on rank ", myRank); flush(io)
  end # of if myRank == 0
  MPI.Barrier(MPI.COMM_WORLD)
  MPI.Finalize()
  return nothing
end


end # of module

## old code
# function MyJulia2_main(InFile1::String, InFile2::String, InFile3::String, InFile4::String,
# 	          TimeLimitInSeconds::Real, ScoringMethod::Int, NetworkModel::String, 
#                   tmstart::UInt64; sol1Dir::String=".")
#     println("MyJulia2 called with:\n",
# 	    "\tInFile1 .......... : ", InFile1, "\n",
# 	    "\tInFile2 .......... : ", InFile2, "\n",
# 	    "\tInFile3 .......... : ", InFile3, "\n",
# 	    "\tInFile4 .......... : ", InFile4, "\n",
# 	    "\tTimeLimitInSeconds : ", TimeLimitInSeconds, "\n",
# 	    "\tScoringMethod .... : ", ScoringMethod, "\n",
# 	    "\tNetworkModel ..... : ", NetworkModel)

#     io = stdout # set up log file
#     OutDir = "./"
    
#     # create MPI manager and workers
#     numprocs = haskey(ENV,"SLURM_NTASKS") ? parse(Int, ENV["SLURM_NTASKS"]) : 1
#         # note: this assumes we use Slurm
#     println(io, "Launching ", numprocs, " contingency solvers.")
#     manager = MPIManager(np=numprocs)
#     addprocs(manager)
#     @eval @everywhere begin
#         push!(LOAD_PATH, dirname(@__FILE__))
#         using SCACOPF, GoUtils, CompetitionInterface, Distributed
# 	println("Solver ", Distributed.myid(), " running at ", Distributed.gethostname(), ".")
#     end
#     println(io, numprocs, " contingency solvers launched and idle.")
    
#     # read instance
#     println(io, "Reading and parsing instance data ...")
#     clock = @elapsed begin
# 	MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches,
# 	switchedshunts, generatordsp, activedsptables, costcurves, governorresponse,
# 	contingencies = readinstance(InFile3, InFile4, InFile2, InFile1)
# 	N, L, T, SSh, G, K, P = 
#             GOfmt2params(MVAbase, buses, loads, fixedbusshunts, generators,
# 			 ntbranches, tbranches, switchedshunts, generatordsp, activedsptables, costcurves,
# 			 governorresponse, contingencies)
#     end
#     println(io, "Data prepared for problem formulation in ", round(clock, digits=2), " secs.")
#     flush(io)

#     @printf("Generators %d   Lines %d   Buses %d  Contingencies %d\n", size(G,1), size(L,1), size(N,1), size(K,1))

#     # compute index and sets for performant formulation
#     L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn, K_outidx = indexsets(N, L, T, SSh, G, K)
    
#     # read base case solution
#     println(io, "Reading base case solution from 'solution1.txt' ...")
#     clock = @elapsed begin
# 	v_n, theta_n, b_s, p_g, q_g = ReadBaseSolution(sol1Dir, MVAbase, N, SSh, G)
# 	_, prodcost, penaltybase = BaseCaseValue(v_n, theta_n, b_s, p_g, q_g,
# 					      N, L, T, SSh, G, P)
# 	SCACOPFSol = SCACOPF.SCACOPFState(SCACOPF.BaseState(v_n, theta_n, b_s,
# 							    p_g, q_g, prodcost,
# 							    penaltybase))
#     end
#     println(io, "Base case solution read in ", round(clock, digits=2), " secs."); 
    
#     # solve contingencies given the base solution in SCACOPFSol.base_case
#     allconttime = @elapsed begin
	
# 	# define options and other auxiliary vars
# 	o = SCACOPF.GoOptions()
# 	o.tmstart = tmstart
# 	o.CouplingMode = :complementarityapprox
# 	o.SmoothingParamCoupling=1e-4
# 	o.ipopt.mu_init = 1
# 	o.ipopt.linear_solver = "ma57"
# 	o.ipopt.tol           = 1e-8
# 	o.ipopt.mu_target     = 1e-9
# 	o.ipopt.print_level = 0
# 	SolverCon = SCACOPF.with_ipopt_optimizer(o)
	
# 	# evaluation function
# 	function EvalContingency(k::Int, StartingPoint=nothing)
# 		println("Solving contingency ", k, " out of ", size(K, 1), ".")
# 		contime = @elapsed begin
# 			ContSol, NewStartingPoint = SCACOPF.SolveContingencyTwoStep(o, N, L,
# 				T, SSh, G, (K[:ConType][k], K[:IDout][k], K[:Contingency][k]), P,
# 				SCACOPFSol, SolverCon, StartingPoint=StartingPoint,
# 				IndexSets=(L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn,
# 				Lin, Tidxn, Tin, SShn, Gn, K_outidx[k]))
# 		@assert ContSol.cont_id == K[:Contingency][k]
# 		end
# 		println("\nContingency ", k, " out of ",  size(K, 1), " solved in ",
# 			round(contime, digits=3), " secs.\n")
# 		if StartingPoint == nothing
# 			return k, ContSol, contime, NewStartingPoint
# 		else
# 			return k, ContSol, contime
# 		end
# 	end
        
# 	# parallel loop 1: execute one instance per worker to initial point
#     	println(io, "First parallel loop -- generating initial contingency solutions:");
#         o.SmoothingParamCoupling=1e-2
# 	contsolsinit = pmap(EvalContingency, CachingPool(workers()), 1:min(numprocs, size(K, 1)))
	
# 	# parallel loop
# 	if size(K, 1) > numprocs
#                 o.SmoothingParamCoupling=1e-4
# 		StartingPoint = contsolsinit[1][4]
# 		println(io, "\n\nParallel contingency evaluation loop:");
# 		contsolsvec = pmap(k -> EvalContingency(k, StartingPoint), CachingPool(workers()), 1:size(K, 1))
# 	end

# 	# map solution vector into ordered dictionary (this only creates pointers)
# 	partime = 0.0
        
# 	ContSolutions = OrderedDict{Int, SCACOPF.ContingencyState}()
# 	for i = 1:length(contsolsinit)
# 		k = contsolsinit[i][1]
# 		ContSol = contsolsinit[i][2]
# 		soltime = contsolsinit[i][3]
#                 ContSolutions[K[:Contingency][k]] = ContSol
#                 partime += soltime
# 	end
# 	if size(K, 1) > numprocs
# 		for i = 1:length(contsolsvec)
# 			k = contsolsvec[i][1]
# 			ContSol = contsolsvec[i][2]
# 			soltime = contsolsvec[i][3]
# 			ContSolutions[K[:Contingency][k]] = ContSol
#                 	partime += soltime
# 		end
# 	end
        
#     end #allconttime
#     println(io, "All contingencies solved in ", round(allconttime, digits=2),
# 		" secs. Parallel efficiency: ", round(100*partime/(allconttime*numprocs), digits=1),
# 		"%." ); flush(io)
	
#     # write contingency solution ("solution2.txt", using GO competition format)
#     println(io, "Writting contingencies solutions ...")
#     writetime = @elapsed begin
#         #! for now just copy sol to the arrays so that the IO code is not touched
#         #! we will revisit this later

#         @assert(length(ContSolutions)==size(K,1))
        
# 	v_nk     = Array{Float64, 2}(undef, size(N, 1), size(K, 1))
# 	theta_nk = Array{Float64, 2}(undef, size(N, 1), size(K, 1))
# 	b_sk     = Array{Float64, 2}(undef, size(SSh, 1), size(K, 1))
# 	p_gk     = Array{Float64, 2}(undef, size(G, 1), size(K, 1))
# 	q_gk     = Array{Float64, 2}(undef, size(G, 1), size(K, 1))
# 	delta_k  = Vector{Float64}(undef, size(K, 1))
        
#         for k=1:size(K,1)
#             ContSol = ContSolutions[K[:Contingency][k]]
#             @assert K[:Contingency][k]==ContSol.cont_id
            
#             v_nk[:,k]     = ContSol.v_n
#             theta_nk[:,k] = ContSol.theta_n
#             b_sk[:,k]     = ContSol.b_s
#             p_gk[:,k]     = ContSol.p_g
#             q_gk[:,k]     = ContSol.q_g
#             delta_k[k]    = ContSol.delta
#         end

# 	writesolution(OutDir,
#                       SCACOPFSol.base_state.v_n, 
#                       SCACOPFSol.base_state.theta_n,
#                       SCACOPFSol.base_state.b_s,
#                       SCACOPFSol.base_state.p_g,
#                       SCACOPFSol.base_state.q_g,  
#                       MVAbase, N, SSh, G, generators, true, 
#                       v_nk, theta_nk, b_sk,
# 		      p_gk, q_gk, delta_k, K, contingencies)
#     end
#     println(io, "Contingencies solutions written in ", round(writetime, digits=2), " secs."); flush(io)
	
#     # remove MPI workers
#     rmprocs(manager)

#     return nothing
# end


# ## master call function for SCACOPF
# function mastercall(o_master::GoOptions, 
#                     InFile1::String, InFile2::String, InFile3::String, InFile4::String,
# 	            TimeLimitInSeconds::Real=3600, ScoringMethod::Union{Nothing, Int}=nothing,
# 	            NetworkModel::Union{Nothing, String}=nothing; 
# 	            OutDir::String="./", LogIO::Union{String, IO}=stdout, NumCont::Union{Nothing, Int}=nothing,
# 	            Seed::Union{Nothing, Int}=nothing, SolverSCACOPF=nothing, SolverCon=SolverSCACOPF)::Nothing
	
#     # set up log file
#     if typeof(LogIO) <: String
# 	io = open(LogIO, "w")
#     else
# 	io = LogIO
#     end
    
#     # read instance
#     println(io, "Reading and parsing instance data ...")
#     clock = @elapsed begin
# 	MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches,
# 	switchedshunts, generatordsp, activedsptables, costcurves, governorresponse,
# 	contingencies = readinstance(InFile3, InFile4, InFile2, InFile1)
# 	N, L, T, SSh, G, K, P = 
#             GOfmt2params(MVAbase, buses, loads, fixedbusshunts, generators,
# 			 ntbranches, tbranches, switchedshunts, generatordsp, activedsptables, costcurves,
# 			 governorresponse, contingencies)
#     end
#     println(io, "Data prepared for problem formulation in ", round(clock, digits=2), " secs.")
#     flush(io)

#     @printf("Generators %d   Lines %d   Buses %d  Contingencies %d\n", size(G,1), size(L,1), size(N,1), size(K,1))

#     if size(L,1)>1999 || size(N,1)>1999
#         NumCont = (ScoringMethod==1 || ScoringMethod==3) ? 1 : 6;
#     end

#     # compute index and sets for performant formulation
#     L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn, K_outidx =
#         indexsets(N, L, T, SSh, G, K)
    
#     # select up to NumCont contingencies in conidx
#     NumCont = min(NumCont, size(K,1))    
#     conidx = selectContingencies(NumCont, G, K)

#     SCACOPFSol = nothing
# 	# solve Challenge 1 problem
# 	println(io, "Solving Challenge 1 problem for the selected instance ...")

#         print(io, "Contingencies considered in SCACOPF: ", K[conidx,:], "\n")

#         println(io, " - Pass 1")
#         o1 = o_master
# 	clock1 = @elapsed begin       
#             ## set desired options
#             o1.CouplingMode = :complementarityapprox#:ignore#:fischerburm
#             o1.SmoothingParamCoupling=1e-3
#             # ipopt
#             o1.ipopt.mu_init = 1
#             o1.ipopt.linear_solver = "ma57"
#             o1.ipopt.tol           = 1e-7
#             o1.ipopt.mu_target     = 5e-9
#             o1.ipopt.max_cpu_time  = TimeLimitInSeconds-3
#             #o1.ipopt.print_user_options = "yes"
            
#             SCACOPFSolver = with_ipopt_optimizer(o1)

#             try
# 	        SCACOPFSol1 = SolveSCACOPF(o1, N, L, T, SSh, G, K[conidx,:], P, 
#                                            SCACOPFSolver,
# 			                   IndexSets=nothing)
#                 SCACOPFSol = SCACOPFSol1
#             catch ex
#                 println(ex)
#                 println("exception in pass 1")
#             end
# 	end
#         println(io, "Pass 1 problem solved in ", round(clock1, digits=2), " secs.")


#         if false
#         println(io, " - Pass 2")
#         o2 = o_master
#         if nothing == SCACOPFSol
#             ## pass1 failed
#             conidx=[] #solve with no conting
#             o2.CouplingMode = :ignore#:complementarityapprox #:fischerburm
#             o2.ipopt.linear_solver = "ma57"
#             o2.ipopt.tol           = 1e-8
#             o2.ipopt.mu_target     = 0
#             o2.ipopt.max_cpu_time  = TimeLimitInSeconds-clock1-8

#         else
#             # pass 1 was succesfully
#             ## set desired options
#             o2.CouplingMode = :complementarityapprox #:fischerburm
#             o2.SmoothingParamCoupling=1e-2
#             # ipopt
#             o2.ipopt.mu_init = 1e-3
#             o2.ipopt.linear_solver = "ma57"
#             o2.ipopt.tol           = 1e-8
#             o2.ipopt.mu_target     = 1e-9
#             o2.ipopt.max_cpu_time  = TimeLimitInSeconds-clock1-8
#         end
            
# 	clock2 = @elapsed begin        
#             SCACOPFSolver = with_ipopt_optimizer(o2)
#             try
#                 SCACOPFSol2 = SolveSCACOPF(o2, N, L, T, SSh, G, K[conidx,:], P, 
#                                            SCACOPFSolver,
# 			                   IndexSets=nothing, StartingPoint=SCACOPFSol)
#                 SCACOPFSol = SCACOPFSol2
                
#             catch ex
#                 println(ex)
#                 println("exception in the second pass")
#             end
# 	end 
#         println(io, "Pass 2 problem solved in ", round(clock2, digits=2), " secs."); flush(io)
#         println(io, "SCACOPF problem solved in a total of ", 
#                 round(clock1+clock2, digits=2), " secs."); flush(io)
#         end #if true or false

# 	# write base case solution ("solution1.txt", using GO competition format)
# 	println(io, "Writting base case solution ...")
# 	clock = @elapsed begin
# 		mkpath(OutDir)
#                 writesolution(OutDir, 
#                               SCACOPFSol.base_state.v_n, 
#                               SCACOPFSol.base_state.theta_n,
#                               SCACOPFSol.base_state.b_s,
#                               SCACOPFSol.base_state.p_g,
#                               SCACOPFSol.base_state.q_g, 
#                               MVAbase, N, SSh, G, generators)
                
# 	end
# 	println(io, "Base case solution written in ", round(clock, digits=2), " secs.")
# 	flush(io)

# 	# check whether we want to post-solve contingencies at this point
# 	if !o_master.PostSolveContingencies
# 		println(io, "\nVerifying solution value ...")
# 		clock = @elapsed begin
# 		    obj, prodcost, penaltycost = BaseCaseValue(SCACOPFSol.base_state.v_n, 
#                                                                SCACOPFSol.base_state.theta_n,
#                                                                SCACOPFSol.base_state.b_s,
#                                                                SCACOPFSol.base_state.p_g,
#                                                                SCACOPFSol.base_state.q_g, 
#                                                                N, L, T, SSh, G, P)
# 		end
# 		println(io, "Solution value recomputed in ", round(clock, digits=2), " secs.")
# 		println(io, "Objective value: \$", round(obj, digits=1))
# 		println(io, "Production cost: \$", round(prodcost, digits=1))
# 		println(io, "Penalty cost:    \$", round(penaltycost, digits=1), "\n")
# 		flush(io)
# 		return nothing
# 	end
#         #
# 	# solve contingencies given the base solution
#         #
# 	println(io, "Solving contingencies ...")
#         ContSolutions = OrderedDict{Int, SCACOPF.ContingencyState}()

# 	allconttime = @elapsed begin

#             #compute index sets for K
#             L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn, K_outidx = indexsets(N, L, T, SSh, G, K)
#             ContSolApprox = nothing
#             cont_list = []
#             append!(cont_list, conidx)
#             append!(cont_list, setdiff(1:size(K,1), conidx))
#             cont_count = 0
#             for k in cont_list
#                 o_cont = o_master
#                 o_cont.CouplingMode = :complementarityapprox #:fischerburm
#                 o_cont.SmoothingParamCoupling=1e-4
#                 # ipopt
#                 o_cont.ipopt.mu_init = 1
#                 o_cont.ipopt.linear_solver = "ma57"
#                 o_cont.ipopt.tol           = 1e-8
#                 o_cont.ipopt.mu_target     = 1e-9
#                 o_cont.ipopt.print_level = 5
                
#                 SolverCon = with_ipopt_optimizer(o_cont)
#                 cont_count = cont_count +1
#                 print(io, "Solving contingency ", cont_count, " with ID=", K[:Contingency][k], " out of ", size(K, 1), "\n")

#                 contime = @elapsed begin

#                     @assert !haskey(ContSolutions,K[:Contingency][k])
# 		    (ContSol, ContSolApprox) = SolveContingencyTwoStep(o_cont, N, L, T, SSh, G, 
#                                                       (K[:ConType][k], K[:IDout][k], K[:Contingency][k]), P,
# 				                      SCACOPFSol, SolverCon, StartingPoint=ContSolApprox,
#                                                       IndexSets=(L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn,
# 				                                 Lin, Tidxn, Tin, SShn, Gn, K_outidx[k]))

#                     @assert ContSol.cont_id == K[:Contingency][k]
#                     ContSolutions[K[:Contingency][k]] = ContSol
# 		end
                
#                 print(io, "Contingency ", k, " solved in ", round(contime, digits=3), " secs.\n\n");flush(io)
#             end

# 	end
# 	println(io, "All contingencies solved in ", round(allconttime, digits=2), " secs."); flush(io)
	
# 	# write contingency solution ("solution2.txt", using GO competition format)
# 	println(io, "Writting contingencies solutions ...")
# 	writetime = @elapsed begin
#             #! for now just copy sol to the arrays so that the IO code is not touched
#             #! we will revisit this later

#             #this will crap out in MPI mode
#             @assert(length(ContSolutions)==size(K,1))

# 	    v_nk     = Array{Float64, 2}(undef, size(N, 1), size(K, 1))
# 	    theta_nk = Array{Float64, 2}(undef, size(N, 1), size(K, 1))
# 	    b_sk     = Array{Float64, 2}(undef, size(SSh, 1), size(K, 1))
# 	    p_gk     = Array{Float64, 2}(undef, size(G, 1), size(K, 1))
# 	    q_gk     = Array{Float64, 2}(undef, size(G, 1), size(K, 1))
# 	    delta_k  = Vector{Float64}(undef, size(K, 1))
            
#             for k=1:size(K,1)
#                 ContSol = ContSolutions[K[:Contingency][k]]
#                 @assert K[:Contingency][k]==ContSol.cont_id

#                 v_nk[:,k]     = ContSol.v_n
#                 theta_nk[:,k] = ContSol.theta_n
#                 b_sk[:,k]     = ContSol.b_s
#                 p_gk[:,k]     = ContSol.p_g
#                 q_gk[:,k]     = ContSol.q_g
#                 delta_k[k]    = ContSol.delta
#             end

# 	    writesolution(OutDir,
#                           SCACOPFSol.base_state.v_n, 
#                           SCACOPFSol.base_state.theta_n,
#                           SCACOPFSol.base_state.b_s,
#                           SCACOPFSol.base_state.p_g,
#                           SCACOPFSol.base_state.q_g,  
#                           MVAbase, N, SSh, G, generators, true, 
#                           v_nk, theta_nk, b_sk,
# 			  p_gk, q_gk, delta_k, K, contingencies)
# 	end
# 	println(io, "Contingencies solutions written in ", round(writetime, digits=2), " secs.")
# 	flush(io)
	
# 	# close log file
# 	if typeof(LogIO) <: String
# 		close(io)
# 	end
	
# 	# return
# 	return nothing
	
# end
