## function to read an instance

function readinstance(RAWfname::T, ROPfname::T, INLfname::T, CONfname::T) where {T <: AbstractString}
    MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches, switchedshunts = readRAW(RAWfname)
    generatordsp, activedsptables, costcurves = readROP(ROPfname)
    governorresponse = readINL(INLfname)
    contingencies = readCON(CONfname)
    return MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches,
    switchedshunts, generatordsp, activedsptables, costcurves, governorresponse,
    contingencies
end

## function to strip IDs of extra stuff

clean_id(str::AbstractString) = strip(replace(str, "'" => ""))

## function to detect starts and ends points of sections

RAWCOMMENT = "@!"

function iscomment(l::AbstractString)::Bool
    return startswith(l, RAWCOMMENT)
end

function next_non_comment(f)
    l = readline(f)
    while iscomment(l)
        l = readline(f)
    end
    return l
end

function isheadortail(l::AbstractString)::Bool
	if length(l) == 0
		return false
	elseif l[1] != '0' && l[1] != ' '
		return false
	elseif length(l) == 1 && l[1] == '0'
		return true
	elseif length(l) >= 2 && l[1:2] == "0 "
		return true
	elseif length(l) >= 3 && l[1:3] == " 0 "	# noncompliant but present in some ROP files
		return true
	else
		return false
	end
end

function sections(filename::AbstractString, casedata::Bool=true, returnheaders::Bool=false)
	f = open(filename, "r")
	sectionstarts = Int[]
	sectionends = Int[]
    l = readline(f)
	if casedata
		append!(sectionstarts, [1 + iscomment(l), 4])
		push!(sectionends, 3 + iscomment(l))
	else
		push!(sectionstarts, 1 + iscomment(l))
	end
	if returnheaders
		headers = String[]
		for i = 1:length(sectionstarts)
			push!(headers, "")
		end
	end
	nlines = 1
	while !eof(f)
        l = readline(f)
		nlines += 1
		if length(l) == 0
			continue
		end
		if isheadortail(l)
			push!(sectionends, nlines-1)
			push!(sectionstarts, nlines+1)
			if returnheaders
				push!(headers, l)
			end
		elseif l[1] == 'Q' && length(l) == 1
			push!(sectionends, nlines-1)
		end
	end
	close(f)
	if length(sectionends) < length(sectionstarts)
		push!(sectionends, nlines)
	end
	if returnheaders
		return sectionstarts, sectionends, headers
	else
		return sectionstarts, sectionends
	end
end

## function to get index of next non-empty subsection

simplifystr(str::String) = lowercase(replace(str, " " => ""))

function getnextsectionidx(currentidx::Int, sectionstarts::AbstractVector{Int},
	sectionends::AbstractVector{Int}, headers::Union{AbstractVector{String}, Nothing}=nothing,
	keyword::Union{String, Nothing}=nothing)::Union{Int, Nothing}
	if headers != nothing
		keyword = simplifystr(keyword)
	end
	for idx = (currentidx+1):length(sectionstarts)
		if sectionends[idx] >= sectionstarts[idx] &&
			(headers == nothing || occursin(keyword, simplifystr(headers[idx])))
			return idx
		end
	end
	return nothing
end

## function to read RAW file
# RAW file contains the following information:
# + System MVA base
# + Buses: id, area, voltage magnitude and angle, lb and ub for voltage in normal and emergency condtions
# + Loads: bus, id, status, active and reactive constant power
# + Fixed bus shunts: bus, id, status, active and reactive power at v=1pu
# + Generators: bus, id, real and reactive power output, lb and ub reactive power, status, lb and ub active power
# + Non-transformer branches: from bus, to bus, circuit, resistance, reactance, susceptance, normal and emergency ratings at v=1pu, status
# + Transformer branches (2-windings): 
# + Switched shunts: bus, status, initial susceptance, steps and susceptance per block

function readRAW(filename::AbstractString)

	# find starting and ending point of each data section
	secstarts, secends = sections(filename, true)
	
	# get MVA base
	f = open(filename, "r")
    l = next_non_comment(f)
	close(f)
	MVAbase = parse(Float64, split(l, ',')[2]) 
	
	# read bus data
    print("Reading buses ... ")
	buses = CSV.read(filename, DataFrame,
		header=[:I,:NAME,:BASKV,:IDE,:AREA,:ZONE,:OWNER,:VM,:VA,:NVHI,:NVLO,:EVHI,:EVLO],
		skipto=secstarts[2], limit=secends[2]-secstarts[2]+1, delim=',',
		types=Dict(1=>Int, 4=>Int, 5=>Int, 8=>Float64, 9=>Float64, 10=>Float64, 11=>Float64,
			12=>Float64, 13=>Float64), ntasks=1, comment=RAWCOMMENT, quotechar='\'')
    println("done. Read ", size(buses, 1), " buses.")
	
	# read load data
	if secends[3] >= secstarts[3]
        print("Reading loads ... ")
		loads = CSV.read(filename, DataFrame,
			header=[:I,:ID,:STATUS,:AREA,:ZONE,:PL,:QL,:IP,:IQ,:YP,:YQ,:OWNER,:SCALE,:INTRPT],
			skipto=secstarts[3], limit=secends[3]-secstarts[3]+1, delim=',',
			types=Dict(1=>Int, 2=>String, 3=>Int, 6=>Float64, 7=>Float64), ntasks=1,
            comment=RAWCOMMENT)
		loads[!,:ID] = clean_id.(loads[!,:ID])
        println("done. Read ", size(loads, 1), " buses.")
	else
		loads = DataFrame()
	end
	
	# read fixed bus shunt data
	if secends[4] >= secstarts[4]
        print("Reading shunts ... ")
		fixedbusshunts = CSV.read(filename, DataFrame,
			header=[:I,:ID,:STATUS,:GL,:BL],
			skipto=secstarts[4], limit=secends[4]-secstarts[4]+1, delim=',',
			types=Dict(1=>Int, 2=>String, 3=>Int, 4=>Float64, 5=>Float64), ntasks=1,
            comment=RAWCOMMENT)
		fixedbusshunts[!,:ID] = clean_id.(fixedbusshunts[!,:ID])
        println("done. Read ", size(fixedbusshunts, 1), " shunts.")
	else
		fixedbusshunts = DataFrame()
	end
	
	# generator data
    print("Reading generators ... ")
	generators = CSV.read(filename, DataFrame,
		header=[:I,:ID,:PG,:QG,:QT,:QB,:VS,:IREG,:MBASE,:ZR,:ZX,:RT,:XT,:GTAP,:STAT,:RMPCT,
			:PT,:PB,:O1,:F1,:O2,:F2,:O3,:F3,:O4,:F4,:WMOD,:WPF],
		skipto=secstarts[5], limit=secends[5]-secstarts[5]+1, delim=',',
		types=Dict(1=>Int, 2=>String, 3=>Float64, 4=>Float64, 5=>Float64, 6=>Float64, 15=>Int,
			17=>Float64, 18=>Float64), ntasks=1, comment=RAWCOMMENT)
	generators[!,:ID] = clean_id.(generators[!,:ID])
    println("done. Read ", size(generators, 1), " generators.")

	# non-transformer branch data
	if secends[6] >= secstarts[6]
        print("Reading branches ... ")
		ntbranches = CSV.read(filename, DataFrame,
			header=[:I,:J,:CKT,:R,:X,:B,:RATEA,:RATEB,:RATEC,:GI,:BI,:GJ,:BJ,:ST,:MET,:LEN,
				:O1,:F1,:O2,:F2,:O3,:F3,:O4,:F4],
			skipto=secstarts[6], limit=secends[6]-secstarts[6]+1, delim=',',
			types=Dict(1=>Int, 2=>Int, 3=>String, 4=>Float64, 5=>Float64, 6=>Float64,
				7=>Float64, 9=>Float64, 14=>Int), ntasks=1, comment=RAWCOMMENT)
		ntbranches[!,:CKT] = clean_id.(ntbranches[!,:CKT])
        println("done. Read ", size(ntbranches, 1), " branches.")
	else
		ntbranches = DataFrame()
	end
	
	# transformer data
	if secends[7] >= secstarts[7]
        print("Reading transformers ... ")
		tbranches = readtransformerdata(filename, secstarts[7], secends[7])
        println("done. Read ", size(tbranches, 1), " transformers.")
	else
		tbranches = DataFrame()
	end
	
	# switched shunt data
	if secends[18] >= secstarts[18]
        print("Reading switched shunts ... ")
		switchedshunts = CSV.read(filename, DataFrame,
			header=[:I,:MODSW,:ADJM,:STAT,:VSWHI,:VSWLO,:SWREM,:RMPCT,:RMIDNT,:BINIT,
				:N1,:B1,:N2,:B2,:N3,:B3,:N4,:B4,:N5,:B5,:N6,:B6,:N7,:B7,:N8,:B8],
			skipto=secstarts[18], limit=secends[18]-secstarts[18]+1, delim=',',
			types=Dict(1=>Int, 4=>Int, 10=>Float64, 11=>Float64, 12=>Float64, 13=>Float64,
				14=>Float64, 15=>Float64, 16=>Float64, 17=>Float64, 18=>Float64,
				19=>Float64, 20=>Float64, 21=>Float64, 22=>Float64, 23=>Float64,
				24=>Float64, 25=>Float64, 26=>Float64), ntasks=1, comment=RAWCOMMENT,
                maxwarnings=0)
        println("done. Read ", size(switchedshunts, 1), " switched shunts.")
	else
		switchedshunts = DataFrame()
	end
	
	# return data frames with RAW data
	return MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches, switchedshunts

end

## function to read ROP file
# ROP file contains the following information:
# + Generator dispatch: correspondance between generators and dispatch tables
# + Active power dispatch tables: correspondance between dispatch tables and cost curves
# + Cost curves: piecewise linear cost curves, described as pairs (production_i, cost_i)
#				 polynomial cost curves, described as polynomial coeffient (costquad, 
#				 costlin, cost)

function readROP(filename::AbstractString)
	
    # if file is empty, return here
    if read(filename, String) == ""
        return DataFrame(), DataFrame(), DataFrame()
    end
	
    # find starting and ending point of each data section
	secstarts, secends, headers = sections(filename, false, true)
    
	# generator dispatch data
	idx = getnextsectionidx(1, secstarts, secends, headers, "Generator Dispatch")
	if secends[idx] >= secstarts[idx]
		generatordsp = CSV.read(filename, DataFrame,
			header=[:BUS,:GENID,:DISP,:DSPTBL],
			skipto=secstarts[idx], limit=secends[idx]-secstarts[idx]+1,
			quotechar='\'', types=Dict(1=>Int, 2=>String, 4=>Int), ntasks=1)
		generatordsp[!,:GENID] = clean_id.(generatordsp[!,:GENID])
	else
		generatordsp = DataFrame()
	end
	
	# active power dispatch tables
	idx = getnextsectionidx(idx, secstarts, secends, headers, "Active Power Dispatch")
	if secends[idx] >= secstarts[idx]
		activedsptables = CSV.read(filename, DataFrame,
			header=[:TBL,:PMAX,:PMIN,:FUELCOST,:CTYP,:STATUS,:CTBL],
			skipto=secstarts[idx], limit=secends[idx]-secstarts[idx]+1,
			quotechar='\'', types=Dict(1=>Int, 7=>Int), ntasks=1)
	else
		activedsptables = DataFrame()
	end
	
	if activedsptables[1, :CTYP] == 1
		# polynomial and exponentail cost curves
		idx = getnextsectionidx(idx, secstarts, secends, headers, "Polynomial Cost")
		if secends[idx] >= secstarts[idx]
			costcurves = readpolynomialcostcurves(filename, secstarts[idx], secends[idx])
		else
			costcurves = DataFrame()
		end
	else
		# piecewise linear cost curves
		idx = getnextsectionidx(idx, secstarts, secends, headers, "Piece-wise Linear Cost")
		if secends[idx] >= secstarts[idx]
			costcurves = readcostcurves(filename, secstarts[idx], secends[idx])
		else
			costcurves = DataFrame()
		end
	end
	
	# return data frames with ROP data
	return generatordsp, activedsptables, costcurves 

end

## function to read INL file
# INL file contains a single section describing governoor response (generator contingency re-dispatch factor)

function readINL(filename::AbstractString)::DataFrame
	
	# find starting and ending point of each data section (unique in this case)
	secstarts, secends = sections(filename, false)

	# read and return governor response data
	if secends[1] >= secstarts[1]
		governorresponse = CSV.read(filename, DataFrame,
			header=[:I,:ID,:H,:PMAX,:PMIN,:R,:D],
			skipto=secstarts[1], limit=secends[1]-secstarts[1]+1,
			types=Dict(1=>Int, 2=>String, 6=>Float64), ntasks=1)
		governorresponse[!,:ID] = clean_id.(governorresponse[!,:ID])
	else
		governorresponse = DataFrame()
	end
	return governorresponse
end

## function to read CON file
# CON file contains a single section describing all contingencies that can occur in the system
# file is read line-by-line and the result is returned as a data frame

abstract type Contingency end

struct GeneratorContingency <: Contingency
	Bus::Int
	Unit::String
end

struct TransmissionContingency <: Contingency
	FromBus::Int
	ToBus::Int
	Ckt::String
end

function readCON(filename::AbstractString)::DataFrame
	
	# read contingency data
	f = open(filename, "r")
	labels = String[]
	ctypes = Vector{Vector{Symbol}}()
	cons = Vector{Vector{Contingency}}()
	emptycons = String[]
	con_key::Union{String, Nothing} = nothing
	while !eof(f)
		l = readline(f)
		if l == "END"
			break
		end
		if l[1:11] != "CONTINGENCY"
			error("expected contingency start line, found: ", l)
		end
		conname = split(l)[2]
		info = split(readline(f))
		if info[1] != "REMOVE" && info[1] != "OPEN" && info[1] != "END"
			error("expected REMOVE, OPEN or END, found: ", info[1])
		elseif info[1] == "END"
			push!(emptycons, conname)
			continue
		else
			cons_k = Vector{Contingency}()
			ctype = Vector{Symbol}()
			while true
				if info[1] == "REMOVE"		# generator contingency
					push!(cons_k, GeneratorContingency(parse(Int, info[6]), info[3]))
					push!(ctype, :Generator)
					con_key = info[1]
				elseif info[1] == "OPEN"	# branch contingency
					push!(cons_k, TransmissionContingency(parse(Int, info[5]),
									parse(Int, info[8]), strip(info[10])))
					push!(ctype, :Branch)
					con_key = info[1]
				elseif info[1] == "END"
					push!(labels, conname)
					push!(cons, cons_k)
					push!(ctypes, ctype)
					break
				else
					error("expected REMOVE, OPEN or END, found: ", info[1])
				end
				info = split(readline(f))
			end
		end
	end
	if length(emptycons) > 0
		@warn(string("contingency registers ", emptycons, " are empty and they will be ignored."))
	end

	# put contingency data in data frame and return
	contingencies = DataFrame([labels, ctypes, cons], [:LABEL, :CTYPE, :CON])
	return contingencies
	
end

## function to read transformer data

function readtransformerdata(filename::AbstractString, startline::Int, endline::Int)::DataFrame
	
	# check that lines number make sense
	if mod(endline - startline + 1, 4) != 0
		error("number of lines must be a multiple of 4")
	end
	
	# collect transformers information, row-by-row
	rawinfo = Array{String, 2}(undef, div(endline - startline + 1, 4), 43)
	f = open(filename, "r")
	for j = 1:(startline-1)
		readline(f)
	end
	i = startline
	k = 1
	while i <= endline
		row = String[]
		for j=1:4
			l = readline(f)
			append!(row, strip.(split(l, ',')))
		end
        try
            rawinfo[k,:] = row
        catch e
            @show k
            @show i
            @show row
            throw(e)
        end
		k += 1
		i += 4
	end
	close(f)
	
	# form a data frame with the collected data and return
	transformers = DataFrame(rawinfo, [:I,:J,:K,:CKT,:CW,:CZ,:CM,:MAG1,:MAG2,
		:NMETR,:NAME,:STAT,:O1,:F1,:O2,:F2,:O3,:F3,:O4,:F4,:VECGRP,:R12,:X12,
		:SBASE12,:WINDV1,:NOMV1,:ANG1,:RATA1,:RATB1,:RATC1,:COD1,:CONT1,:RMA1,
		:RMI1,:VMA1,:VMI1,:NTP1,:TAB1,:CR1,:CX1,:CNXA1,:WINDV2,:NOMV2])
	colnames = names(transformers)
	intcols = [1:2;12]
	floatcols = [8:9;22:23;25;27:28;30;42]
	stringcols = [4,11]
	for col in intcols
		transformers[!,colnames[col]] =
			try
				parse.(Int, transformers[!,colnames[col]])
			catch
				error("failed to parse column ", col, " (", colnames[col], ") as Int.")
			end
	end
	for col in floatcols
		transformers[!,colnames[col]] =
			try
				parse.(Float64, transformers[!,colnames[col]])
			catch
				error("failed to parse column ", col, " (", colnames[col], ") as Float64.")
			end
	end
	for col in stringcols
		transformers[!,colnames[col]] = replace.(transformers[!,colnames[col]], Ref("'" => ""))
	end
	transformers[!,:CKT] = clean_id.(transformers[!,:CKT])
	return transformers
	
end

## function to read piecewise linear cost functions

function readcostcurves(filename::AbstractString, startline::Int, endline::Int)::DataFrame

	# collect cost curve information, row-by-row
	lbtl = Int[]
	label = String[]
	npairs = Int[]
	xi = Vector{Vector{Float64}}()
	yi = Vector{Vector{Float64}}()
	f = open(filename, "r")
	for i = 1:(startline-1)
		readline(f)
	end
	i = startline
	while i <= endline
		l = readline(f)
		header = strip.(split(l, ','))
		if length(header) != 3
			error("cost curve should start with a 3 field line, got: ", l)
		end
		push!(lbtl, parse(Int, header[1]))
		push!(label, replace(header[2], "'" => ""))
		push!(npairs, parse(Int, header[3]))
		x = Float64[]
		y = Float64[]
		for j = 1:npairs[end]
			xy = parse.(Float64, strip.(split(readline(f), ',')))
			if j == 1 || xy[1] > x[end]
				push!(x, xy[1])
				push!(y, xy[2])
			end
		end
		push!(xi, x)
		push!(yi, y)
		i += 1 + npairs[end]
	end
	close(f)
	
	# place cost curve information in a data frame and return
	costcurves = DataFrame([lbtl, label, npairs, xi, yi], [:LTBL,:LABEL,:NPAIRS,:Xi,:Yi])
	return costcurves
	
end

## function to read polynomial and exponential cost functions

function readpolynomialcostcurves(filename::AbstractString, startline::Int, endline::Int)::DataFrame

	# collect cost curve information, row-by-row
	pbtl = Int[]
	label = String[]
	cost = Vector{Float64}()
	costlin = Vector{Float64}()
	costquad = Vector{Float64}()
	costexp = Vector{Float64}()
	expn = Vector{Float64}()
	f = open(filename, "r")
	for i = 1:(startline-1)
		readline(f)
	end
	i = startline
	while i <= endline
		l = readline(f)
		header = strip.(split(l, ','))
		if length(header) != 7
			error("cost curve should start with a 7 field line, got: ", l)
		end
		push!(pbtl, parse(Int, header[1]))
		push!(label, replace(header[2], "'" => ""))
		push!(cost, parse(Float64, header[3]))
		push!(costlin, parse(Float64, header[4]))
		push!(costquad, parse(Float64, header[5]))
		push!(costexp, parse(Float64, header[6]))
		push!(expn, parse(Float64, header[7]))

		i += 1
	end
	close(f)
	
	# place cost curve information in a data frame and return
	costcurves = DataFrame([pbtl, label, cost, costlin, costquad, costexp, expn], 
	             [:PLTBL,:LABEL,:COST,:COSTLIN,:COSTQUAD,:COSTEXP,:EXPN])
	return costcurves
	
end

## function to go from bus code to explanatory symbol

function IDE_to_type(type_code::Int)::Symbol
    if type_code == 1
        return :PQ
    elseif type_code == 2
        return :PV
    elseif type_code == 3
        return :SWING
    elseif type_code == 4
        return :DISCONNECTED
    else
        throw(DomainError("type_code must be an integer in {1, 2, 3, 4}"))
    end
end

## function to compute conductance and susceptance

IMPEDANCE_TOL = 1E-4

function conductance_g(r::Float64, x::Float64;
                       atol::Float64=IMPEDANCE_TOL, fallback_inf::Float64=Inf)
    z = r + im * x
    if abs(z) <= atol
        if r > 0
            return Inf
        elseif r < 0
            return -Inf
        else
            @assert r == 0
            return fallback_inf
        end
    else
        return real(inv(z))
    end
end

function susceptance_b(r::Float64, x::Float64;
                       atol::Float64=IMPEDANCE_TOL, fallback_inf::Float64=-Inf)
    z = r + im * x
    if abs(z) <= atol
        if x > 0
            return -Inf
        elseif x < 0
            return Inf
        else
            @assert x == 0
            return fallback_inf
        end
    else
        return imag(inv(z))
    end
end

## function to convert 0 rates to Inf (capacity 0 means Inf in RAW)

if_zero_inf(x::Float64) = x == 0 ? Inf : x

## function to transform GO competition tables into usable tables

function GOfmt2params(MVAbase::Float64, buses::DataFrame, loads::DataFrame,                         # RAW
                      fixedbusshunts::DataFrame, generators::DataFrame, ntbranches::DataFrame,      # RAW
                      tbranches::DataFrame, switchedshunts::DataFrame,                              # RAW
                      generatordsp::Union{DataFrame, Nothing}=nothing,                              # ROP
                      activedsptables::Union{DataFrame, Nothing}=nothing,                           # ROP
                      costcurves::Union{DataFrame, Nothing}=nothing,                                # ROP
	                  governorresponse::Union{DataFrame, Nothing}=nothing,                          # INL
                      contingencies::Union{DataFrame, Nothing}=nothing;                             # CON
                      enforce_bounds_on_x0::Bool=true)

    
    # check input types are coherent
    if typeof(generatordsp) != typeof(activedsptables) ||
       typeof(activedsptables) != typeof(costcurves)
        throw(ArgumentError("ROP params must be all be 'nothing' or all specified in as 'DataFrame'"))
    end
    if typeof(governorresponse) != typeof(contingencies)
        throw(ArgumentError("INL and CON params must be all be 'nothing' or all specified in as 'DataFrame'"))
    end
    
	# buses
	N = DataFrame(Any[buses[!,:I], buses[!,:AREA], buses[!, :BASKV],
		zeros(Float64, size(buses,1)), zeros(Float64, size(buses,1)),
		zeros(Float64, size(buses,1)), zeros(Float64, size(buses,1)),
		buses[!,:NVLO], buses[!,:NVHI], buses[!,:EVLO], buses[!,:EVHI],
		buses[!,:VM], buses[!,:VA]*pi/180, map(IDE_to_type, buses[!,:IDE])],
		[:Bus, :Area, :Vbase, :Pd, :Qd, :Gsh, :Bsh, :Vlb, :Vub, :EVlb, :EVub, :v0, :theta0, :Type])
	if size(loads, 1) > 0
		BusLoad = indexin(loads[!,:I], N[!,:Bus])
		for l = 1:size(loads, 1)
			if BusLoad[l] == nothing
				error("bus ", loads[l,:I], " of load ", l, " not found.")
			end
			if loads[l,:STATUS] == 1
				N[BusLoad[l],:Pd] += loads[l,:PL]/MVAbase
				N[BusLoad[l],:Qd] += loads[l,:QL]/MVAbase
			end
		end
	end
	if size(fixedbusshunts, 1) > 0
		BusShunt = indexin(fixedbusshunts[!,:I], N[!,:Bus])
		for fbsh = 1:size(fixedbusshunts, 1)
			if BusShunt[fbsh] == nothing
				error("bus ", fixedbusshunts[fbsh,:I], " of fixed shunt ", fbsh, " not found.")
			end
			if fixedbusshunts[fbsh,:STATUS] == 1
				N[BusShunt[fbsh],:Gsh] += fixedbusshunts[fbsh,:GL]/MVAbase
				N[BusShunt[fbsh],:Bsh] += fixedbusshunts[fbsh,:BL]/MVAbase
			end
		end
	end
	
	# non-transformer branches
	if size(ntbranches, 1) > 0
	    BusLineFrom = indexin(ntbranches[!,:I], N[!,:Bus])
        BusLineTo = indexin(ntbranches[!,:J], N[!,:Bus])
		activeidxs = findall(i -> ntbranches[i,:ST] != 0 && 
                                  N[BusLineFrom[i],:Type] != :DISCONNECTED &&
                                  N[BusLineTo[i],:Type] != :DISCONNECTED,
                             1:size(ntbranches, 1))
		activelines = view(ntbranches, activeidxs, :)
		L = DataFrame(Any[activeidxs, activelines[!,:I],
			activelines[!,:J], activelines[!,:CKT],
            conductance_g.(activelines[!,:R], activelines[!,:X]),
            susceptance_b.(activelines[!,:R], activelines[!,:X]),
			#activelines[!,:R]./(activelines[!,:R].^2 + activelines[!,:X].^2),
			#-activelines[!,:X]./(activelines[!,:R].^2 + activelines[!,:X].^2),
			activelines[!,:B],
            if_zero_inf.(activelines[!,:RATEA]./MVAbase),
            if_zero_inf.(activelines[!,:RATEC]./MVAbase)],
			[:Line, :From, :To, :CktID, :G, :B, :Bch, :RateBase, :RateEmer])
	else
		L = DataFrame(Any[Int[], Int[], Int[], String[],
			Float64[], Float64[], Float64[], Float64[], Float64[]],
			[:Line, :From, :To, :CktID, :G, :B, :Bch, :RateBase, :RateEmer])
	end
	
	# transformers
	if size(tbranches, 1) > 0
	    BusTrfFrom = indexin(tbranches[!,:I], N[!,:Bus])
        BusTrfTo = indexin(tbranches[!,:J], N[!,:Bus])
		activeidxs = findall(i -> tbranches[i,:STAT] != 0 && 
                                  N[BusTrfFrom[i],:Type] != :DISCONNECTED &&
                                  N[BusTrfTo[i],:Type] != :DISCONNECTED,
                             1:size(tbranches, 1))
		activetrafos = view(tbranches, activeidxs, :)
		T = DataFrame(Any[activeidxs, activetrafos[!,:I], activetrafos[!,:J],
			activetrafos[!,:CKT], activetrafos[!,:MAG1], activetrafos[!,:MAG2],
            conductance_g.(activetrafos[!,:R12], activetrafos[!,:X12]),
            susceptance_b.(activetrafos[!,:R12], activetrafos[!,:X12]),
			#activetrafos[!,:R12]./(activetrafos[!,:R12].^2 + activetrafos[!,:X12].^2),
			#-activetrafos[!,:X12]./(activetrafos[!,:R12].^2 + activetrafos[!,:X12].^2),
			activetrafos[!,:WINDV1]./activetrafos[!,:WINDV2], activetrafos[!,:ANG1]*pi/180,
			if_zero_inf.(activetrafos[!,:RATA1]./MVAbase),
            if_zero_inf.(activetrafos[!,:RATC1]./MVAbase)],
			[:Transformer, :From, :To, :CktID, :Gm, :Bm, :G, :B, :Tau, :Theta, :RateBase, :RateEmer])
	else
		T = DataFrame(Any[Int[], Int[], Int[], String[],
			Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]],
			[:Transformer, :From, :To, :CktID, :Gm, :Bm, :G, :B, :Tau, :Theta, :RateBase, :RateEmer])
	end

	# switched shunts
	if size(switchedshunts, 1) > 0
        BusSShunt = indexin(switchedshunts[!,:I], N[!, :Bus])
        activeidxs = findall(i -> switchedshunts[i,:STAT] != 0 &&
                                  N[BusSShunt[i],:Type] != :DISCONNECTED,
                             1:size(switchedshunts, 1))
		activess = view(switchedshunts, activeidxs, :)
		SSh = DataFrame(Any[activeidxs, activess[!,:I],
			Vector{Float64}(undef, length(activeidxs)),
			Vector{Float64}(undef, length(activeidxs)),
			activess[!,:BINIT]./MVAbase], [:SShunt, :Bus, :Blb, :Bub, :b0])
		for ssh = 1:length(activeidxs)
			Blb = 0.0
			Bub = 0.0
			for i = 1:8
                if ismissing(activess[ssh, 10+2*i])
                    continue
                end
				Bsec = activess[ssh, 9+2*i]*activess[ssh, 10+2*i]/MVAbase
				if Bsec < 0
					Blb += Bsec
				else
					Bub += Bsec
				end
			end
			SSh[ssh,:Blb] = Blb
			SSh[ssh,:Bub] = Bub
		end
	else
		SSh = DataFrame(Any[Int[], Int[], Float64[], Float64[], Float64[]],
			[:SShunt, :Bus, :Blb, :Bub, :b0])
	end
	
	# generators -- RAW
    BusGen = indexin(generators[!,:I], N[!, :Bus])
    activeidxs = findall(i -> generators[i,:STAT] !=0 &&
                              N[BusGen[i],:Type] != :DISCONNECTED,
                         1:size(generators, 1))
	activegens = view(generators, activeidxs, :)
	G = DataFrame(Any[activeidxs, activegens[!,:I], activegens[!,:ID],
		activegens[!,:PB]./MVAbase, activegens[!,:PT]./MVAbase, activegens[!,:QB]./MVAbase,
		activegens[!,:QT]./MVAbase, activegens[!,:PG]./MVAbase, activegens[!,:QG]./MVAbase],
		[:Generator, :Bus, :BusUnitNum, :Plb, :Pub, :Qlb, :Qub, :p0, :q0])
	
	# generators -- ROP
	if isnothing(generatordsp)
        @assert isnothing(activedsptables) && isnothing(costcurves)
        @warn "no information on generator costs provided, assuming cost is 0 for all generators"
		G[!,:COST]     = Vector{Float64}(undef, size(G,1))
		G[!,:COSTLIN]  = Vector{Float64}(undef, size(G,1))
		G[!,:COSTQUAD] = Vector{Float64}(undef, size(G,1))
		for g = 1:size(G, 1)
			G[g,:COST]     = 0.0
			G[g,:COSTLIN]  = 0.0
			G[g,:COSTQUAD] = 0.0
		end
        G[!,:CostPi] = Vector{Vector{Float64}}(undef, size(G,1))
        G[!,:CostCi] = Vector{Vector{Float64}}(undef, size(G,1))
        for g = 1:size(G, 1)
            G[g,:CostPi] = Float64[G[g,:Plb], G[g,:Pub]]
            G[g,:CostCi] = Float64[0.0, 0.0]
        end
		# Structure needs the key CTYP
		G[!, :CTYP] = Vector{Int}(undef, size(G, 1))
		G[!, :CTYP] .= 2
    else
        gdspix = indexin(string.(G[!,:Bus], ":", G[!,:BusUnitNum]),
            string.(generatordsp[!,:BUS], ":", generatordsp[!,:GENID]))
        gdsptbl = generatordsp[!,:DSPTBL][gdspix]
        gctbl = activedsptables[!,:CTBL][indexin(gdsptbl, activedsptables[!,:TBL])]
		G[!, :CTYP] = Vector{Int}(undef, size(G, 1))
		G[!, :CTYP] .= activedsptables[1, :CTYP]
		if G[1, :CTYP] == 1
			gctblix = indexin(gctbl, costcurves[!,:PLTBL])
			if any(gctblix .== nothing)
				error("there seems to be missing cost curves for generators: ",
					findall(x->x!=0, gctblix .== nothing))
			end
			gctblix = convert(Array{Int64}, gctblix)
			G[!,:COST]     = Vector{Float64}(undef, size(G,1))
			G[!,:COSTLIN]  = Vector{Float64}(undef, size(G,1))
			G[!,:COSTQUAD] = Vector{Float64}(undef, size(G,1))
			for g = 1:size(G, 1)
				G[g,:COST]     = costcurves[gctblix[g],:COST]
				G[g,:COSTLIN]  = costcurves[gctblix[g],:COSTLIN]
				G[g,:COSTQUAD] = costcurves[gctblix[g],:COSTQUAD]
			end
		else
			gctblix = indexin(gctbl, costcurves[!,:LTBL])
			if any(gctblix .== nothing)
				error("there seems to be missing cost curves for generators: ",
					findall(x->x!=0, gctblix .== nothing))
			end
			gctblix = convert(Array{Int64}, gctblix)
			G[!,:CostPi] = Vector{Vector{Float64}}(undef, size(G,1))
			G[!,:CostCi] = Vector{Vector{Float64}}(undef, size(G,1))
			for g = 1:size(G, 1)
				G[g,:CostPi] = costcurves[gctblix[g],:Xi]./MVAbase
				G[g,:CostCi] = costcurves[gctblix[g],:Yi]
			end
		end
    end

	# ---- fixing infeasible initial dispatchs ----
	if enforce_bounds_on_x0
        modifiedstartingpoint = Int[]
        for g = 1:size(G, 1)
            if G[g,:p0] < G[g,:Plb] || G[g,:p0] > G[g,:Pub]
                G[g,:p0] = .5*G[g,:Plb] + .5*G[g,:Pub]
                push!(modifiedstartingpoint, g)
            end
            if G[g,:q0] < G[g,:Qlb] || G[g,:q0] > G[g,:Qub]
                G[g,:q0] = .5*G[g,:Qlb] + .5*G[g,:Qub]
                push!(modifiedstartingpoint, g)
            end
        end
        if length(modifiedstartingpoint) > 0
            unique!(modifiedstartingpoint)
            msg = "generators with infeasible starting points: "
            for g = modifiedstartingpoint
                msg *= string(G[g,:BusUnitNum], "/", G[g,:Bus], " ")
            end
            @warn(msg)
        end
        modifiedstartingpoint = nothing
    end
    
	# ---- fix bad bounds in cost functions ----
	if !(isnothing(generatordsp)) 
		if activedsptables[1, :CTYP] != 1
			modifiedcostfunction = Int[]
			for g = 1:size(G, 1)
				xi = G[g,:CostPi]
				yi = G[g,:CostCi]
				n = length(xi)
				if xi[1] > G[g,:Plb]
					yi[1] = yi[1] + (yi[2] - yi[1])/(xi[2] - xi[1])*(G[g,:Plb] - xi[1])
					xi[1] = G[g,:Plb]
					push!(modifiedcostfunction, g)
				end
				if xi[n] < G[g,:Pub]
					yi[n] = yi[n] + (yi[n] - yi[n-1])/(xi[n] - xi[n-1])*(G[g,:Pub] - xi[n])
					xi[n] = G[g,:Pub]
					push!(modifiedcostfunction, g)
				end
			end
			if length(modifiedcostfunction) > 0
				unique!(modifiedcostfunction)
				msg = "generators with inconsistent cost functions: "
				for g = modifiedcostfunction
					msg *= string(G[g,:BusUnitNum], "/", G[g,:Bus], " ")
				end
				@warn(msg)
			end
			modifiedcostfunction = nothing
		else
			modifiedcostfunction = nothing		
		end
	else
		modifiedcostfunction = nothing		
	end
	
	# generators -- INL
	if isnothing(governorresponse)
        @warn "no information on governor response, assuming only swing-bus generators can respond"
	    swing_gens_idx = findall(x -> x == :SWING, N[!,:Type][indexin(G[!,:Bus], N[!,:Bus])])
        G[!,:alpha] = zeros(Float64, size(G, 1))
        G[swing_gens_idx, :alpha] .= 1.0
    else    
		if G[1,:CTYP] != 1   
			if size(governorresponse, 1) == 0
				@warn "no information on governor response, assuming only swing-bus generators can respond"
				swing_gens_idx = findall(x -> x == :SWING, N[!,:Type][indexin(G[!,:Bus], N[!,:Bus])])
				G[!,:alpha] = zeros(Float64, size(G, 1))
				G[swing_gens_idx, :alpha] .= 1.0
			else
				ggovrespix = indexin(string.(G[!,:Bus], ":", G[!,:BusUnitNum]),
					string.(governorresponse[!,:I], ":", governorresponse[!,:ID]))
				if any(ggovrespix .== nothing)
					error("there seems to be missing participation factors for generators: ",
						findall(x->x!=0, ggovrespix .== nothing))
				end
				ggovrespix = convert(Array{Int64}, ggovrespix)
				G[!,:alpha] = governorresponse[ggovrespix,:R]
			end
		end
    end
	
	# contingencies
	if isnothing(contingencies)
        # @warn "no information on contingencies"
        K = DataFrame(Any[Int64[], Symbol[], Int64[], String[]],
						[:Contingency,:ConType,:IDout, :Label])
    else
        K = DataFrame(Any[collect(1:size(contingencies, 1)),
			Vector{Vector{Symbol}}(undef, size(contingencies, 1)),
			Vector{Union{Vector{Int64}, Nothing}}(nothing, size(contingencies, 1)),
			Vector{String}(undef, size(contingencies, 1))],
            [:Contingency,:ConType,:IDout,:Label])
        if size(contingencies, 1) > 0
            missingcon_k = Int64[]
            missingcon_el = Int64[]
			for k = 1:size(contingencies, 1)
				con_types = Symbol[]
				idout = Int64[]
				con_type_k = contingencies[k, :CTYPE]
				con_k = contingencies[k, :CON]
				gencon = findall(con_type_k .== :Generator)
				txcon = findall(con_type_k .== :Branch)
				searchstr = String[]
				for el = 1:length(con_type_k)
					if el in gencon
						push!(searchstr, string.(con_k[el].Bus, ":", con_k[el].Unit))
					else
						push!(searchstr, string.(con_k[el].FromBus, ":", 
												con_k[el].ToBus, ":", con_k[el].Ckt))
					end
				end
				gix = indexin(searchstr, string.(G[!,:Bus], ":", G[!,:BusUnitNum]))
				lix = indexin(searchstr, string.(L[!,:From], ":", L[!,:To], ":", L[!,:CktID]))
				trix = indexin(searchstr, string.(T[!,:From], ":", T[!,:To], ":", T[!,:CktID]))
				for el = 1:length(con_type_k)
					if lix[el] != nothing
						push!(idout, L[lix[el],:Line])
						push!(con_types, :Line)
					elseif trix[el] != nothing
						push!(idout, T[trix[el],:Transformer])
						push!(con_types, :Transformer)
					elseif gix[el] != nothing
						push!(idout, G[gix[el],:Generator])
						push!(con_types, :Generator)
					else
						push!(missingcon_k, k)
						push!(missingcon_el, el)
					end
				end
				K[k,:IDout] = idout
				K[k,:ConType] = con_types
				K[k,:Label] = contingencies[k, :LABEL]
			end
            if length(missingcon_k) > 0
                msg = "found inconsistent contingency registers: "
                missingcontbl = view(contingencies, missingcon_k, :)
				missing_type = Symbol[]
				for i = 1:length(missingcon_el)
					push!(missing_type, missingcontbl[i, :CTYPE][missingcon_el[i]])
				end
                gencon = findall(missing_type .== :Generator)
				searchstr = String[]
                if length(gencon) > 0
                    genoff = view(generators, findall(x->x==0, tbranches[!,:STAT]), :)
					for i = 1:length(gencon)
	                    push!(searchstr, string.(missingcontbl[gencon[i],:CON][missingcon_el[gencon[i]]].Bus, ":", 
                        missingcontbl[gencon[i],:CON][missingcon_el[gencon[i]]].Unit))
					end
                    gix = indexin(searchstr, string.(genoff[!,:I], ":", genoff[!,:ID]))
                    for i = 1:length(gencon)
                        gcon = missingcontbl[gencon[i],:CON][missingcon_el[gencon[i]]]
                        msg *= string(" ", missingcontbl[gencon[i],:LABEL], ":",
                            gcon.Bus, "/", gcon.Unit)
                        if gix[i] != nothing
                            msg *= ":STAT=0"
                        else
                            msg *= ":missing"
                        end
                    end
                end
                txcon = findall(missing_type .== :Branch)
                if length(txcon) > 0
                    linoff = view(ntbranches, findall(x->x==0, ntbranches[!,:ST]), :)
                    troff = view(tbranches, findall(x->x==0, tbranches[!,:STAT]), :)
					for i = 1:length(txcon)
						push!(searchstr, string.(missingcontbl[txcon[i],:CON][missingcon_el[txcon[i]]].FromBus, ":",
											missingcontbl[txcon[i],:CON][missingcon_el[txcon[i]]].ToBus, ":",
											missingcontbl[txcon[i],:CON][missingcon_el[txcon[i]]].Ckt))
					end
                    lix = indexin(searchstr, string.(linoff[!,:I], ":", linoff[!,:J], ":", linoff[!,:CKT]))
                    trix = indexin(searchstr, string.(troff[!,:I], ":", troff[!,:J], ":", troff[!,:CKT]))
                    for i = 1:length(txcon)
                        tcon = missingcontbl[txcon[i],:CON][missingcon_el[txcon[i]]]
                        msg *= string(" ", missingcontbl[txcon[i],:LABEL], ":",
                            tcon.FromBus, "/", tcon.ToBus, "/", tcon.Ckt)
                        if lix[i] != nothing || trix[i] != nothing
                            msg *= ":STAT=0"
                        else
                            msg *= ":missing"
                        end
                    end
                end
                msg *= ". These contingencies will be ignored while solving SCACOPF."
                @warn(msg)
				for i = size(K, 1):-1:1
					if 	length(K[i, :ConType]) == 0
	                	deleteat!(K, i)
					end
				end
            end
			@assert all(K[!,:IDout] .!= nothing)
			K[!,:IDout] = convert(Vector{Vector{Int}}, K[!,:IDout])
        end
    end
	
	# penalties
	P = DataFrame(Any[Symbol[:P,:Q,:S],
		[[2, 50, Inf]./MVAbase, [2, 50, Inf]./MVAbase, [2, 50, Inf]./MVAbase],
		[[1E3, 5E3, 1E6].*MVAbase, [1E3, 5E3, 1E6].*MVAbase, [1E3, 5E3, 1E6].*MVAbase]],
		[:Slack, :Quantities, :Penalties])
	#P = DataFrame(Any[Symbol[:P,:Q,:S],
	#	[[Inf]./MVAbase, [Inf]./MVAbase, [Inf]./MVAbase],
	#	[[1E3].*MVAbase, [1E3].*MVAbase, [1E3].*MVAbase]],
	#	[:Slack, :Quantities, :Penalties])
	
    # filter out disconnected buses before returning
    N = N[N[!,:Type] .!= :DISCONNECTED, :]
	
	# return params for SCACOPF
	return N, L, T, SSh, G, K, P
	
end

## function to read and parse an instance

function instancefilenames(instancedir::String, maxnup::Int=3)
        extensions = ["raw", "rop", "inl", "con"]
        exfiles = String[]
        for ex in extensions
                nup = 0
                exfile = ""
                fdir = instancedir
                while nup < maxnup
                        exfile = joinpath(fdir, "case."*ex)
                        isfile(exfile) && break
                        nup += 1
                        fdir *= "/.."
                end
                nup == maxnup && error("case.", ex, " not found")
                push!(exfiles, exfile)
        end
        return tuple(exfiles...)
end

function ParseInstance(dir::String, maxnup::Int=3)
	return ParseInstance(instancefilenames(dir, maxnup)...)
end

function ParseInstance(rawfile::String, ropfile::String, inlfile::String, confile::String)
	MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches,
		switchedshunts, generatordsp, activedsptables, costcurves, governorresponse,
		contingencies = readinstance(rawfile, ropfile, inlfile, confile)
	return GOfmt2params(MVAbase, buses, loads, fixedbusshunts, generators,
		ntbranches, tbranches, switchedshunts, generatordsp, activedsptables, costcurves,
		governorresponse, contingencies)
end
