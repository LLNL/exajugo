#__precompile__()

module SolutionWriterMP

## external modules used

using SparseArrays, Printf, DataFrames, JuMP, CSV, XLSX, Statistics



## elements to be exported

export writesolutionMP, writesolutionBid, writesolutionPlant, writesolutionSlack, writeloss, writepara, writesolutionAllBuses

## function to write solution block

function writesolutionblock(io::IO, v_n::T, theta_n::T, b_s::T, p_g::T, q_g::T,
	MVAbase::Float64, N::DataFrame, SSh::DataFrame, G::DataFrame,
	generators::DataFrame)::Nothing where {T <: AbstractVector{Float64}}
	
	# compute reactive compensaton at each bus
	SSh_Nidx = convert(Vector{Int}, indexin(SSh[!, :Bus], N[!, :Bus]))
	bcsn = spzeros(size(N, 1))
	for ssh = 1:size(SSh, 1)
		bcsn[SSh_Nidx[ssh]] += b_s[ssh]
	end
	bcsn *= MVAbase
	
	# write bus section
	@printf(io, "--bus section\ni, v(p.u.), theta(deg), bcs(MVAR at v = 1 p.u.)\n")
	for n = 1:size(N, 1)
		@printf(io, "%d, %.20f, %.20f, %.20f\n", N[!, :Bus][n], v_n[n], 180/pi*theta_n[n], bcsn[n])
	end
	
	# write generator section
	gmap = zeros(Int, size(generators, 1))
	for g = 1:size(G, 1)
		gmap[G[!, :Generator][g]] = g
	end
	@printf(io, "--generator section\ni, id, p(MW), q(MW)\n")
	for gi = 1:size(generators, 1)
		if gmap[gi] == 0
			@printf(io, "%d, \'%s\', 0, 0\n", generators[!, :I][gi], generators[!, :ID][gi])
		else
			g = gmap[gi]
			@printf(io, "%d, \'%s\', %.12f, %.12f\n", G[!, :Bus][g], G[!, :BusUnitNum][g],
				MVAbase*p_g[g], MVAbase*q_g[g])
		end
	end
	
	return
	
end

## function to write solution (solution1.txt)

function writesolutionMP(OutDir::String, v_n::U, theta_n::U, b_s::U, p_g::U, q_g::U, MVAbase::Float64, 
	N::T, SSh::T, G::T, generators::DataFrame, timestep::Vector{Int64}) ::Nothing where {T <: Vector{Any}, U <: JuMP.Containers.DenseAxisArray}
	
	# write base case solution
	f = open(OutDir * "/solution1.txt", "w")
    for ts = timestep
        @printf(f,"ts: %d\n",ts)
	    writesolutionblock(f, v_n[ts,:], theta_n[ts,:], b_s[ts,:], p_g[ts,:], q_g[ts,:], MVAbase, N[ts], SSh[ts], G[ts], generators)
    end
	close(f)
	
end

## function to compute reactive compensaton at each bus

function MVA(MVAbase::Float64,N::DataFrame, SSh::DataFrame,b_s::AbstractVector{Float64})
	SSh_Nidx = convert(Vector{Int}, indexin(SSh[!, :Bus], N[!, :Bus]))
	bcsn = spzeros(size(N, 1))
	for ssh = 1:size(SSh, 1)
		bcsn[SSh_Nidx[ssh]] += b_s[ssh]
	end
	bcsn *= MVAbase
	return bcsn
end


## function to to write solution block for all ts and n 

function writesolutionblockBid(io::IO, bcsn::SparseVector{Float64, Int64}, MVAbase::Float64, timestep::Vector{Int64}, N::T, G::T, 
	v_n::U, theta_n::U, p_g::U, q_g::U, generators::DataFrame, plant::DataFrame,
	plantpower::Array{Float64}, LMPs::JuMP.Containers.DenseAxisArray{Float64},
	LMPsN::JuMP.Containers.DenseAxisArray{Float64})::Nothing where {T <: Vector{Any}, U <: JuMP.Containers.DenseAxisArray}

	# write bus section
	@printf(io, "--bus section\ni, t, v(p.u.), theta(deg), bcs(MVAR at v = 1 p.u.)\n")
	for n = 1:size(N[1], 1)
		for ts = timestep
			@printf(io, "%d, %d, %.20f, %.20f, %.20f\n", N[ts][!, :Bus][n], ts, v_n[ts,n], 180/pi*theta_n[ts,n], bcsn[n])
		end
	end
		
	# write generator section
	gmap = zeros(Int, size(generators, 1))
	for g = 1:size(G[1], 1)
		gmap[G[1][!, :Generator][g]] = g
	end
	@printf(io, "--generator section\ni, t, id, p(MW), q(MW)\n")
	for gi = 1:size(generators, 1)
		for ts = timestep
			if gmap[gi] == 0
				@printf(io, "%d, %d, \'%s\', 0, 0\n", generators[!, :I][gi], ts, generators[!, :ID][gi])
			else
				g = gmap[gi]
				@printf(io, "%d, %d, \'%s\', %.12f, %.12f\n", G[ts][!, :Bus][g], ts, G[ts][!, :BusUnitNum][g],
				MVAbase*p_g[ts,g], MVAbase*q_g[ts,g])
			end
		end
	end

	@printf(io, "--section total gen\nt, p(MW)\n")
	for ts = timestep
		pg = 0.0
		for gi = 1:size(generators, 1)
			if gmap[gi] != 0
				g = gmap[gi]
				pg += MVAbase*p_g[ts,g]
			end
		end
		@printf(io,"%d, %.12f\n", ts, pg)
	end
	for ts = timestep[end]
		pmin = 0.0
		pmax = 0.0
		for i = 1:size(G[1],1)
			pmin += MVAbase*G[ts][!,:Plb][i]
			pmax += MVAbase*G[ts][!,:Pub][i]
		end
		@printf(io, "--section gen min max\npmin(MW), pmax(MW)\n")
		@printf(io,"%.12f, %.12f\n", pmin, pmax)
	end

	@printf(io, "--section total demand\nt, p(MW)\n")
	for ts = timestep
		pd = 0.0
		for n = 1:size(N[1], 1)
			pd += MVAbase*N[ts][!, :Pd][n]
		end
		@printf(io,"%d, %.12f\n", ts, pd)
	end

	@printf(io, "--section LMPs\ni, t, LMP(USD/MWh)\n")
	nbusn = findall( x -> !(x in plant[!,:plantid]), N[1][!, :Bus])

	for n = eachindex(nbusn)
		for ts = timestep
			@printf(io, "%d, %d, %.12f\n", N[1][!, :Bus][nbusn[n]], ts, LMPsN[ts,n])
		end
	end

	
end

## function to write MPACOPF bidding solution 

function writesolutionBid(OutDir::String, MVAbase::Float64, N::T, SSh::T, G::T, generators::DataFrame, plant::DataFrame,
	timestep::Vector{Int64}, v_n::U, theta_n::U, b_s::U, p_g::U, q_g::U, prodcost::Float64, basepen::Float64, plantrev::Float64, 
	plantpower::Array{Float64}, LMPs::JuMP.Containers.DenseAxisArray{Float64}, 
	LMPsN::JuMP.Containers.DenseAxisArray{Float64}, DBratio::Float64, RatioBidding::Float64, Day::Int64) ::Nothing where {T <: Vector{Any}, U <: JuMP.Containers.DenseAxisArray}

	
	bcsn = MVA(MVAbase, N[1], SSh[1], b_s[1,:])
	f = open(OutDir * "/solutionBid_$(timestep[end])h_M10OZ$(Day)D$(RatioBidding).txt", "w")
	# f = open(OutDir * "/solutionBid_$(lpad(timestep[end],2,"0")).txt", "w")
	@printf(f, "prodcost(USD), basepen(USD), plantrev(USD), DBratio\n")
	@printf(f, "%.20f, %.20f, %.20f, %.20f\n", prodcost, basepen, plantrev, DBratio)
	writesolutionblockBid(f, bcsn, MVAbase, timestep, N, G, v_n, theta_n, p_g, q_g, generators, plant, plantpower, LMPs, LMPsN)
	close(f)
end

## function to write plant solution to CSV files
function writesolutionPlant(OutDir::String, N::T, timestep::Vector{Int64}, plantpower::Array{Float64}, plant::DataFrame,
	LMPs::JuMP.Containers.DenseAxisArray{Float64}, DBratio::Float64, RatioBidding::Float64, Day::Int64) where {T <: Vector{Any}}
		

	nbus = findall( x -> x in plant[!,:plantid], N[1][!, :Bus])
	# println("nbus ",nbus)

	fpb = open(OutDir * "/solutionPlantbus.txt", "w")
	@printf(fpb, "busID,time,plantpower,LMP\n")
	for n = eachindex(nbus)
		for ts = timestep
			@printf(fpb, "%d, %d, %.20f, %.20f\n", N[ts][!, :Bus][nbus[n]], ts, plantpower[ts,n], LMPs[ts,n])
		end
	end

	close(fpb)

	filename2 = joinpath(OutDir,"solutionPlantbus.txt")
	PlantbusSol = CSV.read(filename2, DataFrame)
	rm(filename2)


	CSV.write(OutDir * "/Plantbus_$(timestep[end])h_M10OZ$(Day)D$(RatioBidding).CSV", PlantbusSol);

end

## function to write plant solution to CSV files
function writesolutionAllBuses(OutDir::String, N::T, timestep::Vector{Int64}, plant::DataFrame, Genbus::U, Llossbus::U, 
	LMPs::U, LMPsN::U, RatioBidding::Float64, Day::Int64) where {T <: Vector{Any}, U <: JuMP.Containers.DenseAxisArray{Float64}}

	nbus = findall( x -> x in plant[!,:plantid], N[1][!, :Bus])
	nbusn = findall( x -> !(x in plant[!,:plantid]), N[1][!, :Bus])

	fpb1 = open(OutDir * "/solutionLMP.txt", "w")
	@printf(fpb1, "busID,time,LMP\n")
	for n = eachindex(nbus)
		for ts = timestep
			@printf(fpb1, "%d, %d, %.20f\n", N[ts][!, :Bus][nbus[n]], ts, LMPs[ts,n])
		end
	end
	for n = eachindex(nbusn)
		for ts = timestep
			@printf(fpb1, "%d, %d, %.20f\n", N[1][!, :Bus][nbusn[n]], ts, LMPsN[ts,n])
		end
	end

	close(fpb1)

	filename1 = joinpath(OutDir,"solutionLMP.txt")
	LMPsAll = CSV.read(filename1, DataFrame)
	# println(size(LMPsAll))
	LMPsAllBusSorted = sort(LMPsAll, [:busID, :time])
	# println(first(LMPsAllBusSorted, 30))
	LMPsAllBus = LMPsAllBusSorted[!,:LMP]
	LMPsAllBus = reshape(LMPsAllBus, (timestep[end],size(N[timestep[1]], 1)))
	# println(LMPsAllBus[:,2000])
	rm(filename1)

	fpb = open(OutDir * "/solutionAllbus.txt", "w")
	@printf(fpb, "busID,time,Gen,Loss,LMP\n")
	for n = 1:size(N[timestep[1]], 1)
		for ts = timestep
			@printf(fpb, "%d, %d, %.20f, %.20f, %.20f\n", N[ts][!, :Bus][n], ts, Genbus[ts,n], Llossbus[ts,n],LMPsAllBus[ts,n])
		end
	end

	close(fpb)

	filename2 = joinpath(OutDir,"solutionAllbus.txt")
	AllbusSol = CSV.read(filename2, DataFrame)
	rm(filename2)


	CSV.write(OutDir * "/GLPAllbus_$(timestep[end])h_M10OZ$(Day)D$(RatioBidding).CSV", AllbusSol);

end

## function to write slack values to CSV files
function writesolutionSlack(OutDir::String, pmslack::T, ppslack::T, qmslack::T, qpslack::T, lslack::Array{Float64,3}, 
	tslack::Array{Float64,3}, N::U, L::U, TS::U, timestep::Vector{Int64}, RatioBidding::Float64) ::Nothing where {T <: Array{Float64}, U <: Vector{Any}}

	f1 = open(OutDir * "/solutionSlackPQ.txt", "w")
	@printf(f1, "busID,time,Pm,Pp,Qm,Qp\n")
	for n = 1:size(N[timestep[1]], 1)
		for ts = timestep
			@printf(f1, "%d, %d, %.20f, %.20f, %.20f, %.20f\n", N[ts][!, :Bus][n], ts, pmslack[ts,n], 
			ppslack[ts,n], qmslack[ts,n], qpslack[ts,n])
		end
	end
	close(f1)
	filename1 = joinpath(OutDir,"solutionSlackPQ.txt")
	PQslack = CSV.read(filename1, DataFrame)
	rm(filename1)

	f2 = open(OutDir * "/solutionSlackL.txt", "w")
	@printf(f2, "Line,time,L1,L2\n")
	for l = 1:size(L[timestep[1]],1)
		for ts = timestep
			@printf(f2, "%d, %d, %.20f, %.20f\n", L[ts][!,:Line][l], ts, lslack[ts,l,1], lslack[ts,l,2])
		end
	end
	close(f2)
	filename2 = joinpath(OutDir,"solutionSlackL.txt")
	Lslack = CSV.read(filename2, DataFrame)
	rm(filename2)
	
	f3 = open(OutDir * "/solutionSlackT.txt", "w")
	@printf(f3, "Transformer,time,T1,T2\n")
	for t = 1:size(TS[timestep[1]],1)
		for ts = timestep
			@printf(f3, "%d, %d, %.20f, %.20f\n", TS[ts][!,:Transformer][t], ts, tslack[ts,t,1], tslack[ts,t,2])
		end
	end
	close(f3)
	filename3 = joinpath(OutDir,"solutionSlackT.txt")
	Tslack = CSV.read(filename3, DataFrame)
	rm(filename3)

	XLSX.writetable(OutDir * "/Slack_$(timestep[end])h_$(RatioBidding).xlsx", overwrite=true,
	PQs = (collect(DataFrames.eachcol(PQslack)), DataFrames.names(PQslack)),
	Ls = (collect(DataFrames.eachcol(Lslack)), DataFrames.names(Lslack)),
	Ts = (collect(DataFrames.eachcol(Tslack)), DataFrames.names(Tslack))
	)

end


## function to output losses line and transformer power and percent loss 
function writeloss(OutDir::String, Lloss::T, Gtot::T, eta::T, pline::U, ptran::U, gencost::Vector{Float64},  Laconst::DataFrame, Taconst::DataFrame, L::R, TS::R,
	 timestep::Vector{Int64}, RatioBidding::Float64, Day::Int64) where {T <: Array{Float64}, U <: JuMP.Containers.DenseAxisArray, R <:Vector{Any}}

	f1 = open(OutDir * "/solutionLoss.txt", "w")
	@printf(f1, "time,LineLoss,GenTot,GenCost,Eta\n")
	for ts = timestep
		@printf(f1, "%d, %.20f, %.20f, %.20f, %.20f\n", ts, Lloss[ts], Gtot[ts], gencost[ts], eta[ts])
	end
	@printf(f1, "%.3f, %.20f, %.20f, %.20f, %.20f\n", 0 , sum(Lloss), sum(Gtot), sum(gencost), sum(Lloss)/sum(Gtot))
	@printf(f1, "%.3f, %.20f, %.20f, %.20f, %.20f\n", RatioBidding, mean(Lloss), mean(Gtot), mean(gencost), mean(Lloss)/mean(Gtot))
	close(f1)

	filename1 = joinpath(OutDir,"solutionLoss.txt")
	LossSol = CSV.read(filename1, DataFrame)
	rm(filename1)

	f2 = open(OutDir * "/solutionSlackL.txt", "w")
	@printf(f2, "Line,time,L1,L2\n")
	for l = 1:size(L[timestep[1]],1)
		for ts = timestep
			@printf(f2, "%d, %d, %.20f, %.20f\n", L[ts][!,:Line][l], ts, pline[ts,l,1], pline[ts,l,2])
		end
	end
	close(f2)
	filename2 = joinpath(OutDir,"solutionSlackL.txt")
	PL = CSV.read(filename2, DataFrame)
	rm(filename2)
	
	f3 = open(OutDir * "/solutionSlackT.txt", "w")
	@printf(f3, "Transformer,time,T1,T2\n")
	for t = 1:size(TS[timestep[1]],1)
		for ts = timestep
			@printf(f3, "%d, %d, %.20f, %.20f\n", TS[ts][!,:Transformer][t], ts, ptran[ts,t,1], ptran[ts,t,2])
		end
	end
	close(f3)
	filename3 = joinpath(OutDir,"solutionSlackT.txt")
	PT = CSV.read(filename3, DataFrame)
	rm(filename3)

	XLSX.writetable(OutDir * "/Loss_$(timestep[end])h_M10OZ$(Day)D$(RatioBidding).xlsx", overwrite=true,
	Loss = (collect(DataFrames.eachcol(LossSol)), DataFrames.names(LossSol)),
	PLine = (collect(DataFrames.eachcol(PL)), DataFrames.names(PL)),
	PTran = (collect(DataFrames.eachcol(PT)), DataFrames.names(PT)),
	Lactive = (collect(DataFrames.eachcol(Laconst)),DataFrames.names(Laconst)),
	Tactive = (collect(DataFrames.eachcol(Taconst)),DataFrames.names(Taconst))
	)

	# CSV.write(OutDir * "/Loss_$(timestep[end])h_M10OZ$(Day)D$(RatioBidding).CSV", LossSol);

end

function writepara(OutDir::String, G::T, N::T, plant::T, Day::Int64) where {T <: DataFrame}
	CSV.write(OutDir * "/GenTT$(Day)D.CSV", G);
	CSV.write(OutDir * "/BusData$(Day)D.CSV", N);
	CSV.write(OutDir * "/plantData$(Day)D.CSV", plant)

end

end
