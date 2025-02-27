using JuMP, Ipopt
const MOI = JuMP.MathOptInterface
const MOIU = MOI.Utilities

function addstupidcon!(m, x, f=nothing)
	try
		@show f(rand())
	catch
		println("this is very weird ...")
	end
	@NLconstraint(m, f(x) >= 2)
end

function solverandomproblem()
	f(x) = x.^2
	fprime(x) = 2*x
	fprimeprime(x::T) where T = T(2)
	m = Model(with_optimizer(Ipopt.Optimizer))
	JuMP.register(m, :f, 1, f, fprime, fprimeprime)
	@variable(m, xyz >= rand())
	addstupidcon!(m, xyz)
	@NLobjective(m, Min, f(xyz))
	JuMP.optimize!(m)
        if JuMP.primal_status(m) != MOI.FEASIBLE_POINT
		error("solver failed to find a feasible solution")
        end
	return(JuMP.value(xyz))
end
