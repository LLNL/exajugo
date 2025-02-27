#__precompile__()

module SolutionWriterMPbus

## external modules used

using SparseArrays, Printf, DataFrames, JuMP, CSV, XLSX, Statistics



## elements to be exported
export writesolutionAllBuses

## function to write solutions
function writesolutionAllBuses(OutDir::String, timestep::Vector{Int64}, plantpower::U, LMPs::U, Gen::Vector{Float64}, GenCost::Vector{Float64}, numPlant::Int64, Day::Int64, RatioBidding::Float64) where {U <: JuMP.Containers.DenseAxisArray{Float64}}

    fpb = open(OutDir * "/solutionBus.txt", "w")
    @printf(fpb, "time,Gen,Gencost,GencostAvg,LMP,plantpower\n")
    for ts = timestep
        @printf(fpb, "%d, %.20f, %.20f, %.20f, %.20f, %.20f\n", ts, Gen[ts], GenCost[ts], GenCost[ts]/Gen[ts], LMPs[ts], plantpower[ts]/numPlant)
    end

    close(fpb)

	filename = joinpath(OutDir,"solutionBus.txt")
	BusSol = CSV.read(filename, DataFrame)
	rm(filename)

    CSV.write(OutDir * "/BusSol_$(timestep[end])h_M10OZ$(Day)D$(RatioBidding).CSV", BusSol);
end

end