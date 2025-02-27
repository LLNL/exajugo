## plot plant power and LMPs

# modules used 
using Plots, CSV, DataFrames

OutDir = ARGS[1] 

# load output data files
timestep = [5,12,24]
plantdata = []
for i = eachindex(timestep)
    plantdatai = CSV.read(OutDir * "/Plantbus_$(lpad(timestep[i],2,"0")).CSV", DataFrame)
    push!(plantdata, plantdatai)
end
i = 3

plt = plot(plantdata[i][!,:time], plantdata[i][!,:plantpower], xlabel="time (hour)", ylabel="power (MW)", label="power", linewidth=2)
plot!(twinx(), plantdata[i][!,:LMP], c=:red, ylabel="LMP (\$/MWh)", label="price", leg=:bottomleft, linewidth=2)
# plot!(plantdata[2][!,:time], plantdata[2][!,:plantpower])
display(plt)
savefig(plt, OutDir * "/plant$(lpad(timestep[i],2,"0"))h.png")



