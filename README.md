# ExaJuGO: Exascale Julia Grid Optimization

Simple Julia scrips for solving AC power flow, AC optimal power flow, and
security-constrained AC optimal power flow.

These scripts are intended for experimentation with different (possibly, new)
methods, formulations, and settings for solving these power system problem.
Their implementation, therefore, intentionally avoids excessive encapsulation,
which makes other packages difficult to modify by non-developers.

## Installation

The implementation requires *Julia* 1+. After clonning, issue (from a terminal):

```
$ julia -e "using Pkg; Pkg.activate(\".\"); Pkg.instantiate()"
```

from the repository directory to install all dependencies (submodules, Julia
packages, and Ipopt binaries).

## Usage

### Power flow

To solve AC power flow, issue:

```
$ julia /path/to/ACPF.jl /path/to/raw/file.raw /path/to/solution/folder
```

For example, from the repository folder, issuing:

```
$ julia ACPF.jl ./examples/500bus/case.raw ./example_acpf_solution
```

would solve the power flow problem for the 500-bus system in the `examples/500bus`
directory.

### Optimal power flow

To solve AC optimal power flow, issue:

```
$ julia /path/to/ACOPF.jl  /path/to/raw/file.raw /path/to/raw/file.rop /path/to/solution/folder
```

For example, from the repository folder, issuing:

```
$ julia ACOPF.jl ./examples/500bus/case.raw  ./examples/500bus/case.rop ./example_acopf_solution
```

would solve the optimal power flow problem for the 500-bus system in the
`examples/500bus` directory.

Alternatively, you may specify just the directory where `case.raw` and 
`case.rop` files are located, i.e.,

```
$ julia ACOPF.jl ./examples/500bus ./example_acopf_solution
```

which, in this case, would produce the same result as the command above.

### Security-constrained optimal power flow

To solve security-constrained AC optimal power flow, issue:

```
$ julia /path/to/ACOPF.jl  /path/to/raw/file.raw /path/to/raw/file.rop /path/to/raw/file.con 
/path/to/solution/folder
```

For example, from the repository folder, issuing:

```
$ julia SCACOPF.jl ./examples/500bus/case.raw  ./examples/500bus/case.rop ./examples/500bus/case.con ./example_scacopf_solution
```

would solve the optimal power flow problem for the 500-bus system in the
`examples/500bus` directory.

Alternatively, you may specify just the directory where `case.raw`, `case.rop` and
`case.con` files are located, i.e.,

```
$ julia SCACOPF.jl ./examples/500bus ./example_scacopf_solution
```

which, in this case, would produce the same result as the command above.

An additional directory argument can be added to the ACOPF and SCACOPF command pair to specify where to store the output information, including the solution, power flow constraints, objective value, power flow data, and slack variables. For example:

```
$ julia SCACOPF.jl ./examples/500bus ./example_scacopf_solution ./example_scacopf_system
```

## Authors

ExaJuGO is written by Ignacio Aravena (aravenasolis1@llnl.gov), Nai-Yuan Chiang
(chiang7@llnl.gov), and Cosmin G. Petra (petra1@llnl.gov) from LLNL.

## License

ExaJuGO is distributed under the terms of both the MIT license. See [LICENSE](LICENSE) for
details. All new contributions must be made under the MIT license.

## Acknowledgments

ExaJuGO has been developed under the financial support of:

* U.S. Department of Energy, Office of Advanced Scientific Computing Research (ASCR):
  Exascale Computing Program (ECP) and Applied Math Program.
* U.S. Department of Energy, Advanced Research Projects Agency-Energy (ARPA‑E).

LLNL-CODE-2002985
