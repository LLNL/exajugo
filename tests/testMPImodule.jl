__precompile__()

module testMPImodule

using MPI, Distributed, JuMP, Ipopt

export parfun

function parfun(nelems::Int)
	
	# add parallel processes
	manager = MPIManager(np=2)
	addprocs(manager)
	pool = CachingPool(workers())
	
	# load module on workers
	@eval @everywhere begin
		push!(LOAD_PATH, dirname(@__FILE__))
		using testMPImodule
	end
	
	# define and communicate common
	xvec = rand(nelems)
	addition = rand()
	
	function dummy(x::T)::T where {T <: Real}
		println("adding ", addition, " to ", x)
		m = Model(with_optimizer(Ipopt.Optimizer, print_level=0))
		@variable(m, y >= x + addition)
		@objective(m, Min, y)
		JuMP.optimize!(m)
		return JuMP.value(y)
	end
	
	# parallel evaluation
	xvecplussomething = pmap(dummy, pool, xvec)
	
	# remove processes
	rmprocs(manager)
	
	# results
	@show xvec
	@show xvecplussomething
	
end

end
