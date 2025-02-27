#__precompile__()

module InstanceReaderMPbus

## external modules used

using DataFrames, CSV, Statistics

## elements to be exported

export readinstance, GOfmt2params, DELTA, ParseInstance, GOfmt2paramsMPbus

## constants

const DELTA = 0.5

## function to read an instance

function readinstance(RAWfname::T, ROPfname::T, INLfname::T, CONfname::T, CTfname::T, plantfname::T) where {T <: AbstractString}
    MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches, switchedshunts = readRAW(RAWfname)
    generatordsp, activedsptables, costcurves = readROP(ROPfname)
    governorresponse = readINL(INLfname)
    contingencies = readCON(CONfname)
	changetable = readCT(CTfname)
	plant = readplant(plantfname)
    return MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches,
    switchedshunts, generatordsp, activedsptables, costcurves, governorresponse,
    contingencies, changetable, plant
end


## function to detect starts and ends points of sections

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
	sectionstarts = Int[]
	sectionends = Int[]
	if casedata
		append!(sectionstarts, [1,4])
		push!(sectionends, 3)
	else
		push!(sectionstarts, 1)
	end
	if returnheaders
		headers = String[]
		for i = 1:length(sectionstarts)
			push!(headers, "")
		end
	end
	f = open(filename, "r")
	nlines = 0
	for l in eachline(f)
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
	MVAbase = parse(Float64, split(readline(f), ',')[2]) 
	close(f)
	
	# read bus data
	buses = CSV.read(filename,DataFrame;
		header=[:I,:NAME,:BASKV,:IDE,:AREA,:ZONE,:OWNER,:VM,:VA,:NVHI,:NVLO,:EVHI,:EVLO],
		skipto=secstarts[2], limit=secends[2]-secstarts[2]+1, quotechar='\'',
		types=Dict(1=>Int, 5=>Int, 8=>Float64, 9=>Float64, 10=>Float64, 11=>Float64,
			12=>Float64, 13=>Float64))
	
	# read load data
	if secends[3] >= secstarts[3]
		loads = CSV.read(filename,DataFrame;
			header=[:I,:ID,:STATUS,:AREA,:ZONE,:PL,:QL,:IP,:IQ,:YP,:YQ,:OWNER,:SCALE,:INTRPT],
			skipto=secstarts[3], limit=secends[3]-secstarts[3]+1, quotechar='\'',
			types=Dict(1=>Int, 2=>String, 3=>Int, 6=>Float64, 7=>Float64))
		loads[!,:ID] = strip.(loads[!,:ID])
	else
		loads = DataFrame()
	end
	# println("loads RAW ", loads)
	# read fixed bus shunt data
	if secends[4] >= secstarts[4]
		fixedbusshunts = CSV.read(filename,DataFrame;
			header=[:I,:ID,:STATUS,:GL,:BL],
			skipto=secstarts[4], limit=secends[4]-secstarts[4]+1, quotechar='\'',
			types=Dict(1=>Int, 2=>String, 3=>Int, 4=>Float64, 5=>Float64))
		fixedbusshunts[!,:ID] = strip.(fixedbusshunts[!,:ID])
	else
		fixedbusshunts = DataFrame()
	end
	
	# generator data
	generators = CSV.read(filename,DataFrame;
		header=[:I,:ID,:PG,:QG,:QT,:QB,:VS,:IREG,:MBASE,:ZR,:ZX,:RT,:XT,:GTAP,:STAT,:RMPCT,
			:PT,:PB,:O1,:F1,:O2,:F2,:O3,:F3,:O4,:F4,:WMOD,:WPF],
		skipto=secstarts[5], limit=secends[5]-secstarts[5]+1, quotechar='\'',
		types=Dict(1=>Int, 2=>String, 3=>Float64, 4=>Float64, 5=>Float64, 6=>Float64, 15=>Int,
			17=>Float64, 18=>Float64))
	generators[!,:ID] = strip.(generators[!,:ID])

	# non-transformer branch data
	if secends[6] >= secstarts[6]
		ntbranches = CSV.read(filename,DataFrame;
			header=[:I,:J,:CKT,:R,:X,:B,:RATEA,:RATEB,:RATEC,:GI,:BI,:GJ,:BJ,:ST,:MET,:LEN,
				:O1,:F1,:O2,:F2,:O3,:F3,:O4,:F4],
			skipto=secstarts[6], limit=secends[6]-secstarts[6]+1, quotechar='\'',
			types=Dict(1=>Int, 2=>Int, 3=>String, 4=>Float64, 5=>Float64, 6=>Float64,
				7=>Float64, 9=>Float64, 14=>Int))
		ntbranches[!,:CKT] = strip.(ntbranches[!,:CKT])
	else
		ntbranches = DataFrame()
	end
	
	# transformer data
	if secends[7] >= secstarts[7]
		tbranches = readtransformerdata(filename, secstarts[7], secends[7])
	else
		tbranches = DataFrame()
	end
	
	# switched shunt data
	if secends[18] >= secstarts[18]
		switchedshunts = CSV.read(filename,DataFrame;
			header=[:I,:MODSW,:ADJM,:STAT,:VSWHI,:VSWLO,:SWREM,:RMPCT,:RMIDNT,:BINIT,
				:N1,:B1,:N2,:B2,:N3,:B3,:N4,:B4,:N5,:B5,:N6,:B6,:N7,:B7,:N8,:B8],
			skipto=secstarts[18], limit=secends[18]-secstarts[18]+1, quotechar='\'',
			types=Dict(1=>Int, 4=>Int, 10=>Float64, 11=>Float64, 12=>Float64, 13=>Float64,
				14=>Float64, 15=>Float64, 16=>Float64, 17=>Float64, 18=>Float64,
				19=>Float64, 20=>Float64, 21=>Float64, 22=>Float64, 23=>Float64,
				24=>Float64, 25=>Float64, 26=>Float64))
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

function readROP(filename::AbstractString)
	
	# find starting and ending point of each data section
	secstarts, secends, headers = sections(filename, false, true)
	
	# generator dispatch data
	idx = getnextsectionidx(1, secstarts, secends, headers, "Generator Dispatch")
	if secends[idx] >= secstarts[idx]
		generatordsp = CSV.read(filename,DataFrame;
			header=[:BUS,:GENID,:DISP,:DSPTBL],
			skipto=secstarts[idx], limit=secends[idx]-secstarts[idx]+1,
			quotechar='\'', types=Dict(1=>Int, 2=>String, 4=>Int))
		generatordsp[!,:GENID] = strip.(generatordsp[!,:GENID])
	else
		generatordsp = DataFrame()
	end
	
	# active power dispatch tables
	idx = getnextsectionidx(idx, secstarts, secends, headers, "Active Power Dispatch")
	if secends[idx] >= secstarts[idx]
		activedsptables = CSV.read(filename,DataFrame;
			header=[:TBL,:PMAX,:PMIN,:FUELCOST,:CTYP,:STATUS,:CTBL],
			skipto=secstarts[idx], limit=secends[idx]-secstarts[idx]+1,
			quotechar='\'', types=Dict(1=>Int, 7=>Int))
	else
		activedsptables = DataFrame()
	end
	
	# piecewise linear cost curves
	idx = getnextsectionidx(idx, secstarts, secends, headers, "Piece-wise Linear Cost")
	if secends[idx] >= secstarts[idx]
		costcurves = readcostcurves(filename, secstarts[idx], secends[idx])
	else
		costcurves = DataFrame()
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
		governorresponse = CSV.read(filename,DataFrame;
			header=[:I,:ID,:H,:PMAX,:PMIN,:R,:D],
			skipto=secstarts[1], limit=secends[1]-secstarts[1]+1,
			types=Dict(1=>Int, 2=>String, 6=>Float64))
		governorresponse[!,:ID] = strip.(governorresponse[!,:ID])
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
	ctypes = Symbol[]
	cons = Contingency[]
	emptycons = String[]
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
		if info[1] == "REMOVE"		# generator contingency
			push!(labels, conname)
			push!(ctypes, :Generator)
			push!(cons, GeneratorContingency(parse(Int, info[6]), info[3]))
		elseif info[1] == "OPEN"	# branch contingency
			push!(labels, conname)
			push!(ctypes, :Branch)
			push!(cons, TransmissionContingency(parse(Int, info[5]),
				parse(Int, info[8]), strip(info[10])))
		elseif info[1] == "END"
			push!(emptycons, conname)
			continue
		else
			error("expected REMOVE, OPEN or END, found: ", info[1])
		end
		l = readline(f)
		if l != "END"
			error("expected contingency end line, found: ", l)
		end
	end
	if length(emptycons) > 0
		@warn(string("contingency registers ", emptycons, " are empty and they will be ignored."))
	end
	
	# put contingency data in data frame and return
	contingencies = DataFrame([labels, ctypes, cons], [:LABEL, :CTYPE, :CON])
	return contingencies
	
end

## function to read chgtab.CSV file
# chagtab.CSV file contains a change table including the following columns:
# label: time step; prob: probability of change, all 0, not used;
# table: table to be modified, 8, CT_TAREALOAD indicates area-wide bus load change
# row: row # to be modified, it's area # for TX2000 and zone # for IL200
# col: col # to be modified, 4, CT_LOAD_ALL_P modify all loads, real only
# chgtype: type of change, 1, CT_REP replaces old value by value in CT_NEWVAL 
# newval: quantity to use for replacement value
# for more details, refer to scenarios.m files and idx_ct.m

function readCT(filename::AbstractString)
	changetable = CSV.read(filename,DataFrame)
	return changetable
end

function readplant(filename::AbstractString)
	plant = CSV.read(filename,DataFrame)
	return plant
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
		rawinfo[k,:] = row
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
		transformers[!, colnames[col]] =
			try
				parse.(Int, transformers[!, colnames[col]])
			catch
				error("failed to parse column ", col, " (", colnames[col], ") as Int.")
			end
	end
	for col in floatcols
		transformers[!, colnames[col]] =
			try
				parse.(Float64, transformers[!, colnames[col]])
			catch
				error("failed to parse column ", col, " (", colnames[col], ") as Float64.")
			end
	end
	for col in stringcols
		transformers[!, colnames[col]] = replace.(transformers[!, colnames[col]], Ref("'" => ""))
	end
	transformers[!,:CKT] = strip.(transformers[!,:CKT])
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

## function to load bidding plant parameters for scaleup 
function biddingloadScale(plant::DataFrame, scalePlant::Float64)

    # scale up the plant
    plant[!,:pplantub] .*= scalePlant
    plant[!,:pplantlb] .*= scalePlant
    plant[!,:plantdemand] .*= scalePlant
    # calculate ramping limits
    plant[!,:plantRRL] .= (plant[!,:pplantub]-plant[!,:pplantlb]).*plant[!,:RRLratio]./plant[!,:planttau]*60
    # println("plant info ", plant)
    return plant
end

## function to load bidding plant parameters for unit change
function biddingloadUnit(plant::DataFrame, MVAbase::Float64)

		# change unit to MW/MVA
		plant[!,:pplantub] ./= MVAbase
		plant[!,:pplantlb] ./= MVAbase

		plant[!,:plantrev] .*= MVAbase
		plant[!,:plantprod] .*= MVAbase
		# # calculate ramping limits
		# plant[!,:plantRRL] .= (plant[!,:pplantub]-plant[!,:pplantlb]).*plant[!,:RRLratio]./plant[!,:planttau]*60
		# # println("plant info ", plant)
		return plant
end

## function to calculate total load demand 

function LoadDemand(buses::DataFrame, loads::DataFrame, MVAbase::Float64)
	# calculate total load demand
	demandtot = 0.0
	if size(loads, 1) > 0
		BusLoad = indexin(loads[!,:I], buses[!,:I])
		for l = 1:size(loads, 1)
			if BusLoad[l] == nothing
				error("bus ", loads[!,:I][l], " of load ", l, " not found.")
			end
			if loads[!,:STATUS][l] == 1
				demandtot += loads[!,:PL][l]
			end
		end
	end
	demandtot /= MVAbase
    return demandtot    
end
## function to calculate plant demand from total load demand 
## to calculate number of plants that partially/fully replcing load
function demandload(buses::DataFrame, loads::DataFrame, plant::DataFrame, RatioBidding::Float64, MVAbase::Float64) 
    demandtot = LoadDemand(buses, loads, MVAbase)
	# calculate plant demand and # of plants needed
	demandplant_target = demandtot*RatioBidding
	plantloadfix = plant[!, :plantdemand][1]*(1+plant[!, :plantflex][1])/plant[!, :plantprod][1]/24
	num_plant = Int(round(demandplant_target/plantloadfix))
	DBratio = plantloadfix*num_plant/demandtot
	println("demandtot, DBratio, num_plant ", demandtot," ", DBratio, " ", num_plant)

	return num_plant, DBratio, demandtot
end

## calculate annual min load (hour 2259)
## to calculate number of plants that partially/fully replcing load
function mindemandload(changetable::DataFrame, plant::DataFrame, RatioBidding::Float64, MVAbase::Float64) 

	# calculate min total load demand
	# minhour = 1467
	minhour = 2259
	minday = changetable[in(minhour).(changetable.label),:]
	demandtot = sum(minday.newval)
	demandtot /= MVAbase
	# calculate plant demand and # of plants needed
	demandplant_target = demandtot*RatioBidding
	plantloadfix = plant[!, :plantdemand][1]*(1+plant[!, :plantflex][1])/plant[!, :plantprod][1]/24
	num_plant = Int(round(demandplant_target/plantloadfix))
	DBratio = plantloadfix*num_plant/demandtot
	println("min demandtot(MW), DBratio, num_plant ", demandtot*MVAbase," ", DBratio, " ", num_plant)

	return num_plant, DBratio, demandtot

end

## function to decide plant load and fixed load MW/MVA
function plantloadbus(changetable::DataFrame, loads::DataFrame, plant::DataFrame, RatioBidding::Float64, MVAbase::Float64)
    # find # of plants to replace loads, record the actually bidding load ratio to total load
    numPlant, DBratio, demandtot = mindemandload(changetable, plant, RatioBidding, MVAbase)
	# numPlant, DBratio, demandtot = demandload(buses, loads, plant, RatioBidding, MVAbase)
    # println("Total demand load unscaled (MW): ", demandtot*MVAbase)
    # calculate total plant load and fixed load
    PlantLoadall = demandtot*DBratio

    return DBratio, PlantLoadall, numPlant, demandtot
end

# ## function to avoid replacing multiple loads on the same bus
# function duploads(loads::DataFrame, loadsReplace::DataFrame, numplant::Int64)
# 	loadsReplace = unique!(loadsReplace,:I)
# 	sizeL = size(loadsReplace,1)
# 	# println("numplant ", numplant, " unique loads ", sizeL)
# 	if sizeL != numplant
# 		numdiff = 2*numplant - sizeL
# 		println("multiple loads on one bus, need ", numdiff, " rows", " for ", numplant, " buses")
# 		loadsReplace = loads[partialsortperm(loads[!,:PL], 1:numdiff, rev=true),:]
# 		loadsReplace = unique!(loadsReplace,:I)
# 	end
# return loadsReplace
# end

# ## function to replace loads by bidding plants
# function plantloads(buses::DataFrame, loads::DataFrame, plant::DataFrame, RatioBidding::Float64, MVAbase::Float64)

# 	# find # of plants to replace loads, record the actually bidding load ratio to total load
# 	numPlant, DBratio = demandload(buses, loads, plant, RatioBidding, MVAbase)
# 	# make a copy of the loads data and add a column of index for tracking
# 	loadsRep = rename(loads)
# 	loadsRep.rownum = eachindex(loads.PL)
# 	# find the loads for partial replacing (largest numPlant loads)
# 	loadsPartial = loadsRep[partialsortperm(loadsRep[!,:PL], 1:numPlant, rev=true),:]
# 	loadsPartial = duploads(loadsRep, loadsPartial, numPlant)
# 	# println("loadsPartial  ", loadsPartial)

# 	# partial replacement: record the bus ID for plant and remaining load
# 	for n = 1:numPlant
# 		if n == 1
# 			# replace the bus ID if it's the first plant
# 			plant[n,:plantid] = loadsPartial[n,:I]
# 		else
# 			# add a new plant and record the bus ID
# 			push!(plant,plant[1,:])
# 			plant[n,:plantid] = loadsPartial[n,:I]
# 		end
		
# 	end

# 	# sort plant by id 
# 	plant_sorted = sort(plant, [:plantid])

# 	return loadsPartial, plant_sorted, plant, DBratio
# end


## function to transform Go tables and change table to useble time series data tables

# function GOfmt2paramsMP(MVAbase::Float64, buses::DataFrame, loads::DataFrame, fixedbusshunts::DataFrame,
# 	generators::DataFrame, ntbranches::DataFrame, tbranches::DataFrame, switchedshunts::DataFrame,
# 	generatordsp::DataFrame, activedsptables::DataFrame, costcurves::DataFrame,
# 	governorresponse::DataFrame, contingencies::DataFrame, changetable::DataFrame)
function GOfmt2paramsMPbus(MVAbase::Float64, buses::DataFrame, loads::DataFrame, fixedbusshunts::DataFrame,
	generators::DataFrame, ntbranches::DataFrame, tbranches::DataFrame, switchedshunts::DataFrame,
	generatordsp::DataFrame, activedsptables::DataFrame, costcurves::DataFrame,
	governorresponse::DataFrame, contingencies::DataFrame,changetable::DataFrame, plant::DataFrame, RatioBidding::Float64, Tend::Int64, Day::Int64)

	# define vectors of DataFrames
	vN = []; 
	vL = [];
	vT = []; 
	vSSh = []; 
	vG = []; 
	vK = []; 
	vP = []; 

    # define vectors for time series of loads
    vLoad = [];

	# scale up the single plant by 10 times
	plantunit = biddingloadUnit(plant, MVAbase)
    # println("plant unit MVA: ", plantunit)
    plantbid = biddingloadScale(plantunit, 10.0)
    # println("plant unit + scale of 10: ", plantbid)
    # get the total plant load and fixed load
    DBratio, PlantLoadall, numPlants, demandtotal = plantloadbus(changetable, loads, plantbid, RatioBidding, MVAbase)
    plantbidnum = biddingloadScale(plantbid, float(numPlants))
	# # multiple plants on different buses
	# loadsPartial, plantRep, plant_unsorted, DBratio = plantloads(buses, loads, plantbid, RatioBidding, MVAbase)
	
	# for each time step, there are different changes for different row value
	# timestep = unique!(changetable.label)
	# smaller time horizon for testing
	t = Tend
	timestep = collect(1:t)
	println("MP ACOPF for $t hrs")
	flush(stdout)
	# get the scale-back ratio for loads
	# scalebackRatio = scaleback(timestep, Day, changetable, demandtotal, MVAbase)
	# scalebackRatio = 1.0
	# apply load change for all time steps
	for tstep in timestep
		# find the time step of the day selected
		tstep_day = (Day-1)*24+tstep
		# find the changetable rows in this time step containing load change info
		change_tstep = changetable[in([tstep_day]).(changetable.label),:]
		# apply changes to loads in this time step
		loadsCT = applyLoad(loads, change_tstep)
        demandtot = LoadDemand(buses, loadsCT, MVAbase)
        fixedLoadtot = demandtot - PlantLoadall
        push!(vLoad, fixedLoadtot)
		# println("loads after changetable ", loadsCT)
		# # subtract plant power from buses with plant after chgtab
		# if size(loadsPartial,1) > 0
		# 	for n = eachindex(plantRep[!,:plantid])
		# 		Lidx = loadsPartial[!,:rownum][n]
		# 		plantloadfix = plant[!, :plantdemand][n]*(1+plant[!, :plantflex][n])/plant[!, :plantprod][n]/24
		# 		loadsCT[!,:PL][Lidx] = loadsCT[!,:PL][Lidx] - plantloadfix*MVAbase
		# 		# println("load besides plant ",loadsCT[!,:I][Lidx], " ", loadsCT[!,:PL][Lidx])
		# 	end
		# end
		# @assert all(loadsCT[!,:PL] .>= 0)
		# # println("loadsCT ", loadsCT)

	    # call the function from one period ACOPF to minimize changes/effort
        N, L, T, SSh, G, K, P = GOfmt2params(MVAbase, buses, loadsCT, fixedbusshunts, generators,
			ntbranches, tbranches, switchedshunts, generatordsp, activedsptables, costcurves,
			governorresponse, contingencies);

		# save the new values to the array 
		push!(vN, N)
		push!(vL, L)
		push!(vT, T)
		push!(vSSh, SSh)
		push!(vG, G)
		push!(vK, K)
		push!(vP, P)
	end
    println("PlantLoadall (MW): ", PlantLoadall*MVAbase)
    println("Fixed loads at various time steps (MW): ", vLoad[1:timestep[end]].*MVAbase)
    println("plant scaled numPlant: ", plantbidnum)

	return vN, vL, vT, vSSh, vG, vK, vP, timestep, vLoad, plantbidnum, PlantLoadall, DBratio, numPlants

end

## function to calculate the proper scale-back ratio for loads
function scaleback(timestep::Vector{Int64}, Day::Int64, changetable::DataFrame, demandtotal::Float64, MVAbase::Float64)

	demandtotal *= MVAbase*0.8
	AllLoads = zeros(size(timestep,1))
	for tstep in timestep
		# find the time step of the day selected
		tstep_day = (Day-1)*24+tstep
		# find the changetable rows in this time step containing load change info
		change_tstep = changetable[in([tstep_day]).(changetable.label),:]
		# println(change_tstep)
		AllLoads[tstep] = sum(change_tstep.newval)
	end
	
	println("All loads: ", AllLoads)
	scalemean = demandtotal/mean(AllLoads)
	scalemax = demandtotal/minimum(AllLoads)
	scalemin = demandtotal/maximum(AllLoads)
	scaleavg = (scalemin+scalemax)/2
	scaleavgm = demandtotal*2/(minimum(AllLoads)+maximum(AllLoads))
	scalebackavg = sort([scalemean, scaleavg, scaleavgm])
	scaleback = mean(scalebackavg[2:3])
	# println("scale mean: ", scalemean, "scale avg: ", scaleavg, " avgm: ", scaleavgm)

	AllLoadsScaled = AllLoads.*scaleback
	println("All loads scaled: ", AllLoadsScaled)
	println("scale back: ", scaleback, " min: ", scalemin, " max: ", scalemax)

	return scaleback

end

## function to apply load changes to load data
function applyLoad(loads::DataFrame, change_tstep::DataFrame)

	# add index to loads 
	loadsCT = rename(loads)
	loadsCT.rownum = eachindex(loads.PL)
	# make changes area by area 
	for areaidx in change_tstep.row
		# find rows in this area, calculate the load total 
		arealoads = loadsCT[in(areaidx).(loadsCT[!,:AREA]),:]
		loadtot = sum(arealoads[!,:PL])
		# find the corresponding newval and calculate the scale up/down ratio
		CTaidx = findall( x -> x == areaidx, change_tstep.row)
		scale = sum(change_tstep.newval[CTaidx]./loadtot)
		# apply load change row by row in the area using scale
		for rowidx in arealoads[!,:rownum]
			loadsCT[rowidx,:PL] = loadsCT[rowidx,:PL].*scale
		end
	end
	# println("load after: ", loadsCT[!,:PL][10:20])
	#get the new loadsCT table without row numbers
	loadsCT = select!(loadsCT,Not(:rownum))
	return loadsCT
end



## function to transform GO competition tables into usable tables

function GOfmt2params(MVAbase::Float64, buses::DataFrame, loads::DataFrame, fixedbusshunts::DataFrame,
	generators::DataFrame, ntbranches::DataFrame, tbranches::DataFrame, switchedshunts::DataFrame,
	generatordsp::DataFrame, activedsptables::DataFrame, costcurves::DataFrame,
	governorresponse::DataFrame, contingencies::DataFrame)
	
	# buses
	N = DataFrame([buses[!,:I], buses[!,:AREA],
		zeros(Float64, size(buses,1)), zeros(Float64, size(buses,1)),
		zeros(Float64, size(buses,1)), zeros(Float64, size(buses,1)),
		buses[!,:NVLO], buses[!,:NVHI], buses[!,:EVLO], buses[!,:EVHI], buses[!,:VM], buses[!,:VA]*pi/180],
		[:Bus, :Area, :Pd, :Qd, :Gsh, :Bsh, :Vlb, :Vub, :EVlb, :EVub, :v0, :theta0])
	if size(loads, 1) > 0
		BusLoad = indexin(loads[!,:I], N[!,:Bus])
		for l = 1:size(loads, 1)
			if BusLoad[l] == nothing
				error("bus ", loads[!,:I][l], " of load ", l, " not found.")
			end
			if loads[!,:STATUS][l] == 1
				N[!,:Pd][BusLoad[l]] += loads[!,:PL][l]/MVAbase
				N[!,:Qd][BusLoad[l]] += loads[!,:QL][l]/MVAbase
			end
		end
	end
	if size(fixedbusshunts, 1) > 0
		BusShunt = indexin(fixedbusshunts[!,:I], N[!,:Bus])
		for fbsh = 1:size(fixedbusshunts, 1)
			if BusShunt[fbsh] == nothing
				error("bus ", fixedbusshunts[!,:I][fbsh], " of fixed shunt ", fbsh, " not found.")
			end
			if fixedbusshunts[!,:STATUS][fbsh] == 1
				N[!,:Gsh][BusShunt[fbsh]] += fixedbusshunts[!,:GL][fbsh]/MVAbase
				N[!,:Bsh][BusShunt[fbsh]] += fixedbusshunts[!,:BL][fbsh]/MVAbase
			end
		end
	end
	
	# non-transformer branches
	if size(ntbranches, 1) > 0
		activeidxs = findall(x->x!=0, ntbranches[!,:ST])
		activelines = view(ntbranches, activeidxs, :)
		L = DataFrame(Any[activeidxs, activelines[!,:I],
			activelines[!,:J], activelines[!,:CKT],
			activelines[!,:R]./(activelines[!,:R].^2 + activelines[!,:X].^2),
			-activelines[!,:X]./(activelines[!,:R].^2 + activelines[!,:X].^2),
			activelines[!,:B], activelines[!,:RATEA]./MVAbase, activelines[!,:RATEC]./MVAbase],
			[:Line, :From, :To, :CktID, :G, :B, :Bch, :RateBase, :RateEmer])
	else
		L = DataFrame(Any[Int[], Int[], Int[], String[],
			Float64[], Float64[], Float64[], Float64[], Float64[]],
			[:Line, :From, :To, :CktID, :G, :B, :Bch, :RateBase, :RateEmer])
	end

	# transformers
	if size(tbranches, 1) > 0
		activeidxs = findall(x->x!=0, tbranches[!,:STAT])
		activetrafos = view(tbranches, activeidxs, :)
		T = DataFrame(Any[activeidxs, activetrafos[!,:I], activetrafos[!,:J],
			activetrafos[!,:CKT], activetrafos[!,:MAG1], activetrafos[!,:MAG2],
			activetrafos[!,:R12]./(activetrafos[!,:R12].^2 + activetrafos[!,:X12].^2),
			-activetrafos[!,:X12]./(activetrafos[!,:R12].^2 + activetrafos[!,:X12].^2),
			activetrafos[!,:WINDV1]./activetrafos[!,:WINDV2], activetrafos[!,:ANG1]*pi/180,
			activetrafos[!,:RATA1]./MVAbase, activetrafos[!,:RATC1]./MVAbase],
			[:Transformer, :From, :To, :CktID, :Gm, :Bm, :G, :B, :Tau, :Theta, :RateBase, :RateEmer])
	else
		T = DataFrame(Any[Int[], Int[], Int[], String[],
			Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]],
			[:Transformer, :From, :To, :CktID, :Gm, :Bm, :G, :B, :Tau, :Theta, :RateBase, :RateEmer])
	end

	# switched shunts
	if size(switchedshunts, 1) > 0
		activeidxs = findall(x->x!=0, switchedshunts[!,:STAT])
		activess = view(switchedshunts, activeidxs, :)
		SSh = DataFrame(Any[activeidxs, activess[!,:I],
			Vector{Float64}(undef, length(activeidxs)),
			Vector{Float64}(undef, length(activeidxs)),
			activess[!,:BINIT]./MVAbase], [:SShunt, :Bus, :Blb, :Bub, :B0])
		for ssh = 1:length(activeidxs)
			Blb = 0.0
			Bub = 0.0
			for i = 1:8
				Bsec = activess[ssh, 9+2*i]*activess[ssh, 10+2*i]/MVAbase
				if Bsec < 0
					Blb += Bsec
				else
					Bub += Bsec
				end
			end
			SSh[!,:Blb][ssh] = Blb
			SSh[!,:Bub][ssh] = Bub
		end
	else
		SSh = DataFrame(Any[Int[], Int[], Float64[], Float64[], Float64[]],
			[:SShunt, :Bus, :Blb, :Bub, :B0])
	end
	
	# generators -- RAW
	activeidxs = findall(x->x!=0, generators[!,:STAT])
	activegens = view(generators, activeidxs, :)
	G = DataFrame(Any[activeidxs, activegens[!,:I], activegens[!,:ID],
		activegens[!,:PB]./MVAbase, activegens[!,:PT]./MVAbase, activegens[!,:QB]./MVAbase,
		activegens[!,:QT]./MVAbase, activegens[!,:PG]./MVAbase, activegens[!,:QG]./MVAbase],
		[:Generator, :Bus, :BusUnitNum, :Plb, :Pub, :Qlb, :Qub, :p0, :q0])
	
	# generators -- ROP
	gdspix = indexin(string.(G[!,:Bus], ":", G[!,:BusUnitNum]),
		string.(generatordsp[!,:BUS], ":", generatordsp[!,:GENID]))
	gdsptbl = generatordsp[!,:DSPTBL][gdspix]
	gctbl = activedsptables[!,:CTBL][indexin(gdsptbl, activedsptables[!,:TBL])]
	gctblix = indexin(gctbl, costcurves[!,:LTBL])
	if any(gctblix .== nothing)
		error("there seems to be missing cost curves for generators: ",
			findall(x->x!=0, gctblix .== nothing))
	end
	gctblix = convert(Array{Int64}, gctblix)
	G[!,:CostPi] = Vector{Vector{Float64}}(undef, size(G,1))
	G[!,:CostCi] = Vector{Vector{Float64}}(undef, size(G,1))
	for g = 1:size(G, 1)
		G[!,:CostPi][g] = costcurves[!,:Xi][gctblix[g]]./MVAbase
		G[!,:CostCi][g] = costcurves[!,:Yi][gctblix[g]]

		# # change all zero-cost generator costs 
		# if G[!,:CostCi][g] == [0.0, 0.0]
		# 	G[!,:CostCi][g] = [0.0, 100.0]
		# end
	end

	# ---- fixing infeasible initial solutions fixing bad bounds in cost functions ----
	modifiedstartingpoint = Int[]
	modifiedcostfunction = Int[]
	for g = 1:size(G, 1)
		if G[!,:p0][g] < G[!,:Plb][g] || G[!,:p0][g] > G[!,:Pub][g]
			G[!,:p0][g] = .5*G[!,:Plb][g] + .5*G[!,:Pub][g]
			push!(modifiedstartingpoint, g)
		end
		if G[!,:q0][g] < G[!,:Qlb][g] || G[!,:q0][g] > G[!,:Qub][g]
			G[!,:q0][g] = .5*G[!,:Qlb][g] + .5*G[!,:Qub][g]
			push!(modifiedstartingpoint, g)
		end
		xi = G[!,:CostPi][g]
		yi = G[!,:CostCi][g]
		n = length(xi)
		if xi[1] > G[!,:Plb][g]
			yi[1] = yi[1] + (yi[2] - yi[1])/(xi[2] - xi[1])*(G[!,:Plb][g] - xi[1])
			xi[1] = G[!,:Plb][g]
			push!(modifiedcostfunction, g)
		end
		if xi[n] < G[!,:Pub][g]
			yi[n] = yi[n] + (yi[n] - yi[n-1])/(xi[n] - xi[n-1])*(G[!,:Pub][g] - xi[n])
			xi[n] = G[!,:Pub][g]
			push!(modifiedcostfunction, g)
		end
	end
	if length(modifiedstartingpoint) > 0
		unique!(modifiedstartingpoint)
		msg = "generators with infeasible starting points: "
		for g = modifiedstartingpoint
		    msg *= string(G[!,:BusUnitNum][g], "/", G[!,:Bus][g], " ")
		end
		@warn(msg)
	end
	if length(modifiedcostfunction) > 0
		unique!(modifiedcostfunction)
		msg = "generators with inconsistent cost functions: "
		for g = modifiedcostfunction
		    msg *= string(G[!,:BusUnitNum][g], "/", G[!,:Bus][g], " ")
		end
		@warn(msg)
	end
	modifiedstartingpoint = nothing
	modifiedcostfunction = nothing
	
	# generators -- INL
	ggovrespix = indexin(string.(G[!,:Bus], ":", G[!,:BusUnitNum]),
		string.(governorresponse[!,:I], ":", governorresponse[!,:ID]))
	if any(ggovrespix .== nothing)
		error("there seems to be missing participation factors for generators: ",
			findall(x->x!=0, ggovrespix .== nothing))
	end
	ggovrespix = convert(Array{Int64}, ggovrespix)
	G[!,:alpha] = governorresponse[!,:R][ggovrespix]
	
	# contingencies
	K = DataFrame(Any[collect(1:size(contingencies, 1)),
		Vector{Symbol}(undef, size(contingencies, 1)),
		Vector{Union{Int64, Nothing}}(nothing, size(contingencies, 1))],
		[:Contingency,:ConType,:IDout])
	if size(contingencies, 1) > 0
		missingcon = Int[]
		gencon = findall(contingencies[!,:CTYPE] .== :Generator)
		searchstr = string.(collect(gcon.Bus for gcon=contingencies[!,:CON][gencon]),
			":", collect(gcon.Unit for gcon=contingencies[!,:CON][gencon]))
		gix = indexin(searchstr, string.(G[!,:Bus], ":", G[!,:BusUnitNum]))
		for i = 1:length(gencon)
			if gix[i] != nothing
				K[!,:ConType][gencon[i]] = :Generator
				K[!,:IDout][gencon[i]] = G[!,:Generator][gix[i]]
			else
				push!(missingcon, gencon[i])
			end
		end
		txcon = findall(contingencies[!,:CTYPE] .== :Branch)
		searchstr = string.(collect(tcon.FromBus for tcon=contingencies[!,:CON][txcon]), ":",
			collect(tcon.ToBus for tcon=contingencies[!,:CON][txcon]), ":",
			collect(tcon.Ckt for tcon=contingencies[!,:CON][txcon]))
		lix = indexin(searchstr, string.(L[!,:From], ":", L[!,:To], ":", L[!,:CktID]))
		trix = indexin(searchstr, string.(T[!,:From], ":", T[!,:To], ":", T[!,:CktID]))
		for i = 1:length(txcon)
			if lix[i] != nothing
				K[!,:ConType][txcon[i]] = :Line
				K[!,:IDout][txcon[i]] = L[!,:Line][lix[i]]
			elseif trix[i] != nothing
				K[!,:ConType][txcon[i]] = :Transformer
				K[!,:IDout][txcon[i]] = T[!,:Transformer][trix[i]]
			else
				push!(missingcon, txcon[i])
			end
		end
		if length(missingcon) > 0
			msg = "found inconsistent contingency registers: "
			missingcontbl = view(contingencies, missingcon, :)
			gencon = findall(missingcontbl[!,:CTYPE] .== :Generator)
			if length(gencon) > 0
				genoff = view(generators, findall(x->x==0, tbranches[!,:STAT]), :)
				searchstr = string.(collect(gcon.Bus for gcon=missingcontbl[!,:CON][gencon]),
					":", collect(gcon.Unit for gcon=missingcontbl[!,:CON][gencon]))
				gix = indexin(searchstr, string.(genoff[!,:I], ":", genoff[!,:ID]))
				for i = 1:length(gencon)
					gcon = missingcontbl[!,:CON][gencon[i]]
					msg *= string(" ", missingcontbl[!,:LABEL][gencon[i]], ":",
						gcon.Bus, "/", gcon.Unit)
					if gix[i] != nothing
						msg *= ":STAT=0"
					else
						msg *= ":missing"
					end
				end
			end
			txcon = findall(missingcontbl[!,:CTYPE] .== :Branch)
			if length(txcon) > 0
				linoff = view(ntbranches, findall(x->x==0, ntbranches[!,:ST]), :)
				troff = view(tbranches, findall(x->x==0, tbranches[!,:STAT]), :)
				searchstr = string.(
					collect(tcon.FromBus for tcon=missingcontbl[!,:CON][txcon]), ":",
					collect(tcon.ToBus for tcon=missingcontbl[!,:CON][txcon]), ":",
					collect(tcon.Ckt for tcon=missingcontbl[!,:CON][txcon]))
				lix = indexin(searchstr, string.(linoff[!,:I], ":", linoff[!,:J], ":", linoff[!,:CKT]))
				trix = indexin(searchstr, string.(troff[!,:I], ":", troff[!,:J], ":", troff[!,:CKT]))
				for i = 1:length(txcon)
					tcon = missingcontbl[!,:CON][txcon[i]]
					msg *= string(" ", missingcontbl[!,:LABEL][txcon[i]], ":",
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
			deleterows!(K, missingcon)
		end
		@assert all(K[!,:IDout] .!= nothing)
		K[!,:IDout] = convert(Vector{Int}, K[!,:IDout])
	end
	
	# penalties
	P = DataFrame(Any[Symbol[:P,:Q,:S],
		[[2, 50, Inf]./MVAbase, [2, 50, Inf]./MVAbase, [2, 50, Inf]./MVAbase],
		[[1E3, 5E3, 1E6].*MVAbase, [1E3, 5E3, 1E6].*MVAbase, [1E3, 5E3, 1E6].*MVAbase]],
		[:Slack, :Quantities, :Penalties])
	
	# return params for SCACOPF
	return N, L, T, SSh, G, K, P
	
end

## function to read and parse an instance

function instancefilenames(instancedir::String, maxnup::Int=3)
        extensions = ["raw", "rop", "inl", "con", "CSV"]
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

end
