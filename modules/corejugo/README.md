# CoreJuGO: Core Julia Grid Optimization modules

* `instance_reader.jl`: Parser for SC-ACOPF instances in the format of the ARPA-E Grid Optimization
  Competition Challenge 1 (see [here](https://gocompetition.energy.gov/challenges/challenge-1/input-data-format))
* `solution_writer.jl`: Exporter for SC-ACOPF solutions in the format of the ARPA-E Grid Optimization
  Competition Challenge 1 (see [here](https://gocompetition.energy.gov/challenge-1-output-files-and-format))
* `solution_evaluator.jl`: Functions to compute complete solutions starting from state variables
  `v_n`, `theta_n`, `b_s`, `p_g`, `q_g`
* `go_structs.jl`: Structures to hold power grid data and results of SCACOPF subproblems.
