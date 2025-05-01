#__precompile__()

module CoreJuGO

using Printf, Random, DataFrames, CSV, Graphs, JuMP, Ipopt

export SCACOPFdata, enforce_bounds!, GenericContingency, isequal_struct,
       SubproblemSolution, BasecaseSolution, ContingencySolution, SCACOPFsolution,
       add_contingency_solution!,
       get_full_initial_solution, get_full_solution,
       write_solution, write_aux_solution,
       read_base_solution,
       number_of_connected_subsystems, split_on_connected_subsystems

include("go_structs.jl")
include("instance_reader.jl")
include("solution_evaluator.jl")
include("solution_writer.jl")
include("solution_reader.jl")
include("network_graph_analysis.jl")

end
