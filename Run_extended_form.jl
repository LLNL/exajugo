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
    con_filename = "./examples/California/case_Eaton_SD_2.con"
    
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
                                    "linear_solver" => "ma27",
                                    "sb" => "yes")

    # Run SC-ACOPF
    EF_sol = solve_SC_ACOPF(psd, opt, output_dir = mod_output_dir);

    println("\nDone.")

end

# Run the function
run_SCACOPF()