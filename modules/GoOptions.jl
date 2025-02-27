# options struct
# part of module SCACOPF

export GoOptions, with_ipopt_optimizer, ipopt_opt, ipopt_QN, ipopt_optPrint

using JuMP, Ipopt

# these are ipopt optimizer options; field names/types map into ipopt options names/types
mutable struct IpoptGoOptions
    linear_solver::String 
    # ma57_pre_alloc::Float64 #default to 1.05, use 2
    # ma57_automatic_scaling::String #default to no, use yes
    ma97_solve_blas3::String #default to no
    ma97_scaling::String #default to dynamic
    mu_target::Float64
    tol::Float64
    mu_init::Float64 #Ipopt defaults to 0.1, we use 1

    print_level::Int
    max_cpu_time::Float64
    bound_push::Float64
    slack_bound_push::Float64

    mu_linear_decrease_factor::Float64
    mu_superlinear_decrease_power::Float64
    hessian_approximation::String
    dual_inf_tol::Float64

    constr_viol_tol::Float64 #Ipopt defaults to 1e-4, we use 1e-8
    acceptable_tol::Float64
    acceptable_iter::Int #Ipopt defaults to 15, we use 10
    acceptable_constr_viol_tol::Float64 #Ipopt defaults to 0.01, we use 1e-6

    print_user_options::String
    max_iter::Int
    print_frequency_iter::Int
    print_timing_statistics::String

    bound_relax_factor::Float64
end

#this is our options struct; includes optimization solver(s) options
mutable struct GoOptions
    # :approx, :complementarity, :complementarityapprox, :ignore, :fischerburm, 
    # :fischerburmlifted, or :fischerburmsquare
    CouplingMode::Symbol 
    # SmoothApproximations.DEFAULTMU for :approx 
    # or the relaxation parameter for the others
    SmoothingParamCoupling::Float64 
    ThermalLimits::Symbol # :quadr, :cone
    PostSolveContingencies::Bool
    PostSolveType::Symbol # :optimized or :fast
    PenaltyType::Symbol   # :pcwslin or :quadr

    tmstart::UInt64 #absolute time the code was called

    ipopt::IpoptGoOptions
end

GoOptions() = GoOptions(:fischerburm,
                        1e-2,
                        :quadr,
                        false,
                        :optimized,
			:quadr,
                        time_ns(),
                        IpoptGoOptions("ma97", "no", "dynamic", 0., 1e-8, 1, 
                                       5, 24000, 1e-8, 1e-8,
                                       0.4, 1.2, "exact", 1.,
                                       1e-8, 1e-6, 10, 1e-6,
                                       "no", 1000, 1, "no",
                                       0.))

# IPOPT options for solving subproblems
function ipopt_opt()
    ipopt_optimizer = optimizer_with_attributes(
        Ipopt.Optimizer,
        "linear_solver" => "ma57",
        # "mu_init" => Float64(1),
        # "constr_viol_tol" => 1e-8,
        # "acceptable_iter" => 10,
        # "acceptable_constr_viol_tol" => 1e-6,
        # "ma57_pre_alloc" => 2.0,
        # "ma57_automatic_scaling" => "yes",
        "print_level" => 2,
        # "mu_linear_decrease_factor" => 0.4,
        # "mu_superlinear_decrease_power" => 1.25,
        # "fixed_variable_treatment" => "relax_bounds",
        # "bound_relax_factor" => 1e-8,
        "max_iter" => 300,
        "tol" => 1e-6,
        # "hessian_approximation" => "limited-memory",
        # "derivative_test" => "second-order",
        # "nlp_scaling_method" => "none",
        )
    return ipopt_optimizer
end

# IPOPT options for solving MPACOPF, print iterations
function ipopt_optPrint()
    ipopt_optimizer = optimizer_with_attributes(
        Ipopt.Optimizer,
        "linear_solver" => "ma57",
        "mu_init" => Float64(1),
        "constr_viol_tol" => 1e-8,
        "acceptable_iter" => 10,
        "acceptable_constr_viol_tol" => 1e-6,
        "ma57_pre_alloc" => 2.0,
        "ma57_automatic_scaling" => "yes",
        # "print_level" => 2,
        # "mu_linear_decrease_factor" => 0.4,
        # "mu_superlinear_decrease_power" => 1.25,
        # "fixed_variable_treatment" => "relax_bounds",
        # "bound_relax_factor" => 1e-8,
        "max_iter" => 500,
        # "derivative_test" => "first-order",
        # "tol" => 1e-6,
        # "hessian_approximation" => "limited-memory",
        # "derivative_test" => "second-order",
        # "nlp_scaling_method" => "none",
        )
    return ipopt_optimizer
end

# IPOPT options for solving master problem with registered fnc and gradient
function ipopt_QN()
    ipopt_optimizer = optimizer_with_attributes(
        Ipopt.Optimizer,
        "linear_solver" => "ma57",
        "mu_init" => Float64(0.1),
        # "constr_viol_tol" => 1e-8,
        # "acceptable_iter" => 10,
        # "acceptable_constr_viol_tol" => 1e-6,
        "ma57_pre_alloc" => 2.0,
        "ma57_automatic_scaling" => "yes",
        "max_iter" => 50,
        "hessian_approximation" => "limited-memory",
        "mu_strategy" => "monotone",
        "tol" => 1e-6,
        "print_user_options" => "yes",
        # "print_level" => 5,
        # "derivative_test" => "first-order",
        # "derivative_test_tol" => 1e-3,
        # "derivative_test_perturbation" => 1e-7,
        # "nlp_scaling_method" => "none",
        )
    return ipopt_optimizer
end

function with_ipopt_optimizer(o::GoOptions)

    kwargs = Dict{Symbol, Union{Int, Float64, String}}()

    for option in fieldnames(typeof(o.ipopt))
         kwargs[option] =  getfield(o.ipopt, option)
    end
    return with_optimizer(Ipopt.Optimizer; kwargs...)
end
