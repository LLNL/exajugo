

include("hiop.jl")
include("exajugo_call.jl")


function test_TSI(ptrTSIGPmodel, case_path, pf_limit_file)

    TSIGPModel = ptrTSIGPmodel[]

    GPmodel, data, TSI = TSIGPModel.GPmodel, TSIGPModel.data, TSIGPModel.TSI

    m, psd, st_args = TSACOPF(case_path, "solution", pf_limit_file, GPmodel);

    tsicon = TSIConstraint(psd, GPmodel, st_args)
    tsicon_prime = TSIConstraintPrime(psd, GPmodel, st_args)
    tsicon_prime_prime = TSIConstraintPrimePrime(psd, GPmodel, st_args) 
    register(m, :tsicon, 1, (pg,qg) -> tsicon(pg,qg), (pg,qg) -> tsicon_prime(pg,qg), (pg,qg) -> tsicon_prime_prime(pg,qg)) 
    @constraint(m, tsicon( m[:p_g], m[:q_g]) <= 0 )

    JuMP.optimize!(m)

end
