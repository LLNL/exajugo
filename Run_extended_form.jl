using Pkg;
if dirname(PROGRAM_FILE) == ""
    Pkg.activate(".")
	push!(LOAD_PATH, "./modules")
else
    Pkg.activate(dirname(PROGRAM_FILE))
	push!(LOAD_PATH, string(dirname(PROGRAM_FILE), "/modules"))
end

using Ipopt, JuMP, Printf
using SCACOPFSubproblems

function run_SCACOPF()

    # Indicate the files
    raw_filename = "./examples/California/case.raw"
    rop_filename = "./examples/California/case.rop"
    con_filename = "./examples/California/case.con"

    # Eaton wildfire with safe distance (SD) 2 miles with every percentile
    con_filename = "./examples/California/case_Prob_Monument_SD_12.con"
    
    println("Reading raw file "*raw_filename*" , rop file "*rop_filename*", con file "*con_filename*"")
    psd = SCACOPFdata(raw_filename=raw_filename, rop_filename=rop_filename, 
                        con_filename = con_filename)

    # Directory to output all the data 
    output_dir = "./info_dir"

    # Directory to output all of the data for SCACOPF
    scacopf_solution_dir = "./example_scacopf_solution"
    if !ispath(scacopf_solution_dir)
        mkpath(scacopf_solution_dir)
    end

    # The optimizer. Comment out line 38 if HSL ma27 is not provided
    opt = optimizer_with_attributes(Ipopt.Optimizer,
                                    "linear_solver" => "ma57",
                                    "sb" => "yes")

    opt_acopf = optimizer_with_attributes(Ipopt.Optimizer,
                            "linear_solver" => "ma57",
                            "sb" => "yes",
                            "tol" =>  1e-6,
                            "mu_superlinear_decrease_power" =>  1.25,
                            "mu_linear_decrease_factor" =>  0.4,
                            "max_iter" =>  500,
                            "print_user_options"  =>  "yes"
                            )

    opt_con = optimizer_with_attributes(Ipopt.Optimizer,
                            "linear_solver" => "ma57",
                            "sb" => "yes",
                            "tol" =>  1e-6,
                            "mu_superlinear_decrease_power" =>  1.25,
                            "mu_linear_decrease_factor" =>  0.4,
                            "bound_relax_factor" =>  1e-6,
                            "max_iter" =>  500,
                            "fixed_variable_treatment" =>  "relax_bounds",
                            "jacobian_regularization_value" => 1e-10,
                            "inf_pr_output" =>  "internal",
                            "acceptable_dual_inf_tol" =>  0.01,
                            "acceptable_constr_viol_tol" => 1e-6,
                            "acceptable_compl_inf_tol" => 0.01,
                            "acceptable_iter" => 1,
                            "print_user_options"  =>  "yes"
                            )

    # Run basecase
    basecase_solution = solve_basecase(psd, opt_acopf, use_opt = true);
    realized_cost = basecase_solution.base_cost
    
    # Run Contingency subproblem
    for l = 1:length(psd.K.Contingency)
        contingency_solution = solve_contingency(psd, l, basecase_solution, opt_con, minutes_since_base = 15.0, use_huber_like_penalty = false, use_opt = true);
        realized_cost += (.5 / size(psd.K, 1)) * contingency_solution.cont_cost
    end

end

# Run the function
run_SCACOPF()