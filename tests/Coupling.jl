# read command line parameters

NUMTRIALS = parse(Int, ARGS[1])

# load modules

println("\nLoading modules ...")
clock = @elapsed begin
	# external modules
	using Ipopt, JuMP
	const MOI = JuMP.MathOptInterface
	const MOIU = MOI.Utilities
	# internal modules
	push!(LOAD_PATH, joinpath(dirname(PROGRAM_FILE), "../modules"))
	using SCACOPF, SmoothApproximations
end
println("All modules loaded in ", round(clock, digits=2), " secs.")
flush(stdout)

println("\nDefining dependencies ...")
clock = @elapsed begin

# function to test AGC coupling

sclamp(x) = smoothclamp(x)

function testACGcoupling(Plb::Real, Pub::Real, alpha::Real, pvalue::Real,
	deltavalue::Real, CouplingMode::Symbol)
	if CouplingMode == :approx
		referencevalue = Plb+(Pub-Plb)*sclamp((pvalue+alpha*deltavalue-Plb)/(Pub-Plb))
	else
		referencevalue = pvalue+alpha*deltavalue
		if referencevalue < Plb
			referencevalue = Plb
		elseif referencevalue > Pub
			referencevalue = Pub
		end
	end
	m = Model(with_optimizer(Ipopt.Optimizer))
	@variable(m, delta, start=deltavalue)
	@variable(m, Plb <= p <= Pub, start=pvalue)
	@variable(m, Plb <= pk <= Pub, start=referencevalue)
	if CouplingMode == :approx
		JuMP.register(m, :sclamp, 1, sclamp,
			x -> SmoothApproximations.smoothclampprime(x),
			x -> SmoothApproximations.smoothclampprimeprime(x))
	end
	SCACOPF.addAGCconstraint!(m, p, delta, pk, Plb, Pub, alpha, CouplingMode,
		sclamp=sclamp, SetStart=true)
	@objective(m, Min, (p-pvalue)^2 + (delta-deltavalue)^2)
	JuMP.optimize!(m)
        solvestatusprimal = JuMP.primal_status(m)
        solvestatus = JuMP.termination_status(m)
	if (solvestatusprimal != MOI.FEASIBLE_POINT) && (solvestatusprimal != MOI.NEARLY_FEASIBLE_POINT)
		return :IpoptFailed, :Undetermined
	end
	if CouplingMode == :approx
		referencevalue = Plb+(Pub-Plb)*sclamp((JuMP.value(p)+alpha*JuMP.value(delta)-Plb)/(Pub-Plb))
	else
		referencevalue = JuMP.value(p)+alpha*JuMP.value(delta)
		if referencevalue < Plb
			referencevalue = Plb
		elseif referencevalue > Pub
			referencevalue = Pub
		end
	end
	@show CouplingMode, abs(JuMP.value(pk) - referencevalue)/(Pub-Plb)
	@assert abs(JuMP.value(pk) - referencevalue)/(Pub-Plb) < 1E-4 string("failed AGC test with Plb=",
		Plb, ", Pub=", Pub, ", alpha=", alpha, ", p=", pvalue, ", delta=", deltavalue,
		", coupling=", CouplingMode)
	if JuMP.value((p-pvalue)^2 + (delta-deltavalue)^2) >= 1E-6
		return :IpoptSucceded, :Suboptimal
	else
		return :IpoptSucceded, :Optimal
	end
end

# function to test PVPQ coupling

sstep(x) = smoothstep(x)

function testPVPQcoupling(Qlb::Real, Qub::Real, vvalue::Real, vkvalue::Real,
	qkvalue::Real, CouplingMode::Symbol)
	m = Model(with_optimizer(Ipopt.Optimizer))
	@variable(m, v, start=vvalue)
	@variable(m, vk, start=vkvalue)
	@variable(m, Qlb <= qk <= Qub, start=qkvalue)
	if CouplingMode == :approx
		JuMP.register(m, :sstep, 1, sstep,
			x -> SmoothApproximations.smoothstepprime(x),
			x -> SmoothApproximations.smoothstepprimeprime(x))
	end
	SCACOPF.addPVPQconstraint!(m, qk, v, vk, Qlb, Qub, CouplingMode, sstep=sstep, SetStart=true)
	@objective(m, Min, ((v - vk)-(vvalue - vkvalue))^2 + (qk - qkvalue)^2/(Qub-Qlb)^2)
	JuMP.optimize!(m)
        solvestatusprimal = JuMP.primal_status(m)
        solvestatus = JuMP.termination_status(m)
	if (solvestatusprimal != MOI.FEASIBLE_POINT) && (solvestatusprimal != MOI.NEARLY_FEASIBLE_POINT)
		return :IpoptFailed, :Undetermined
	end
	if CouplingMode == :approx
		referencevalue = Qub + (Qlb - Qub)*sstep(JuMP.value(vk - v))
	else
		if JuMP.value(vk-v) > 1E8
			referencevalue = Qlb
		elseif JuMP.value(vk-v) < -1E-8
			referencevalue = Qub
		else
			referencevalue = JuMP.value(qk)
		end
	end
	@show JuMP.value(vk-v), vkvalue-vvalue
	@show JuMP.value(qk), referencevalue
	@show CouplingMode, abs(JuMP.value(qk) - referencevalue)/(Qub-Qlb)
	@assert abs(JuMP.value(qk) - referencevalue)/(Qub-Qlb) < 1E-4 string(
		"failed PV/PQ switching testing with Qlb=", Qlb, ", Qub=", Qub,
		", vvalue=", vvalue, ", vkvalue=", vkvalue, ", qkvalue=", qkvalue,
		", coupling=", CouplingMode)
	if JuMP.value(((v - vk)-(vvalue - vkvalue))^2 + (qk - qkvalue)^2/(Qub-Qlb)^2) >= 1E-6
		return :IpoptSucceded, :Suboptimal
	else
		return :IpoptSucceded, :Optimal
	end
end


end
println("All dependencies defined in ", round(clock, digits=2), " secs.")
flush(stdout)

## run AGC test over randomly generated points

println("\nTesting AGC coupling constraints ...")
clock = @elapsed begin

ipoptfailed = zeros(Int, 4)
ipoptsubopt = zeros(Int, 4)

for i = 1:NUMTRIALS
	Plb = -1000*rand()
	Pub = 1000*rand()
	alpha = 1000*rand()
	p = Plb+(Pub-Plb)*rand()
	delta = 6*rand() - 3
	res = testACGcoupling(Plb, Pub, alpha, p, delta, :approx)
	if res[1] == :IpoptFailed
		ipoptfailed[1] += 1
	end
	if res[2] == :Suboptimal
		ipoptsubopt[1] += 1
	end
	res = testACGcoupling(Plb, Pub, alpha, p, delta, :complementarity)
	if res[1] == :IpoptFailed
		ipoptfailed[2] += 1
	end
	if res[2] == :Suboptimal
		ipoptsubopt[2] += 1
	end
	res = testACGcoupling(Plb, Pub, alpha, p, delta, :fischerburm)
	if res[1] == :IpoptFailed
		ipoptfailed[3] += 1
	end
	if res[2] == :Suboptimal
		ipoptsubopt[3] += 1
	end
	res = testACGcoupling(Plb, Pub, alpha, p, delta, :fischerburmsquare)
	if res[1] == :IpoptFailed
		ipoptfailed[4] += 1
	end
	if res[2] == :Suboptimal
		ipoptsubopt[4] += 1
	end
end

println("\nIpopt failed instances:\n",
	"\tSmoothened      : ", ipoptfailed[1], "\n",
	"\tComplementarity : ", ipoptfailed[2], "\n",
	"\tFischer-Burm    : ", ipoptfailed[3], "\n",
	"\tFischer-Burm^2  : ", ipoptfailed[4])
println("Ipopt suboptimal (local optimal termination) instances:\n",
	"\tSmoothened      : ", ipoptsubopt[1], "\n",
	"\tComplementarity : ", ipoptsubopt[2], "\n",
	"\tFischer-Burm    : ", ipoptsubopt[3], "\n",
	"\tFischer-Burm^2  : ", ipoptsubopt[4])

end
println("\nAGC coupling test finished in ", round(clock, digits=2), " secs.")
flush(stdout)

## run PV/PQ switching test over randomly generated points

println("\nTesting PV/PQ coupling constraints ...")
clock = @elapsed begin

ipoptfailed = zeros(Int, 4)
ipoptsubopt = zeros(Int, 4)

for i = 1:NUMTRIALS
	Qlb = -1000*rand()
	Qub = 1000*rand()
	v = .9+.2*rand()
	qk = Qlb+(Qub-Qlb)*rand()
	if rand() < 0.5

		vk = v
	else
		vk = .5+rand()
		if vk > v
			qk = Qlb
		elseif vk < v
			qk = Qub
		end
	end
	res = testPVPQcoupling(Qlb, Qub, v, vk, qk, :approx)
	if res[1] == :IpoptFailed
		ipoptfailed[1] += 1
	end
	if res[2] == :Suboptimal
		ipoptsubopt[1] += 1
	end
	res = testPVPQcoupling(Qlb, Qub, v, vk, qk, :complementarity)
	if res[1] == :IpoptFailed
		ipoptfailed[2] += 1
	end
	if res[2] == :Suboptimal
		ipoptsubopt[2] += 1
	end
	res = testPVPQcoupling(Qlb, Qub, v, vk, qk, :fischerburm)
	if res[1] == :IpoptFailed
		ipoptfailed[3] += 1
	end
	if res[2] == :Suboptimal
		ipoptsubopt[3] += 1
	end
	res = testPVPQcoupling(Qlb, Qub, v, vk, qk, :fischerburmsquare)
	if res[1] == :IpoptFailed
		ipoptfailed[4] += 1
	end
	if res[2] == :Suboptimal
		ipoptsubopt[4] += 1
	end
end

println("\nIpopt failed instances:\n",
	"\tSmoothened      : ", ipoptfailed[1], "\n",
	"\tComplementarity : ", ipoptfailed[2], "\n",
	"\tFischer-Burm    : ", ipoptfailed[3], "\n",
	"\tFischer-Burm^2  : ", ipoptfailed[4])
println("Ipopt suboptimal (local optimal termination) instances:\n",
	"\tSmoothened      : ", ipoptsubopt[1], "\n",
	"\tComplementarity : ", ipoptsubopt[2], "\n",
	"\tFischer-Burm    : ", ipoptsubopt[3], "\n",
	"\tFischer-Burm^2  : ", ipoptsubopt[4])

end
println("\nPV/PQ coupling test finished in ", round(clock, digits=2), " secs.")
flush(stdout)
