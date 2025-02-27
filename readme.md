# Solve MP AC OPF with Demand Bidding
Implementation of NLP approach for solving multi-period ACOPF with bidding plants

## Software prerequisites
Julia v1.9 is used for this model. All Julia packages can be installed with Pkg.add(). The model uses IPOPT ma57 solver which requires Coin-HSL license. To install the HSL library properly, follow the instructions in https://github.com/coin-or-tools/ThirdParty-HSL and rename/export path as suggested in IPOPT HSL section https://juliapackages.com/p/ipopt. 

## Code organization 
All the updated codes are in the ``implementation`` folder

``TX20000`` folder contains all the inputs needed for case Texas2000

``tests`` folder contains the driver file and an empty con file

``modules`` folder contains all the modules required for excuting

To execute, issue: 
``julia $PATHTO/driverMPACOPF.jl $INSTDIR $OUTDIR $RatioBidding $Timestep $MethodSolve``

where ``$PATHTO`` is the relative path to driverMPACOPF.jl; ``$INSTDIR`` is the input directory containing all the input files; ``$OUTDIR`` is the output directory; ``$RatioBidding`` is the ratio of bidding plants in the system; ``$Timestep`` is the number of time periods; ``$MethodSolve`` is the method used to solve the problem.

Currently, ``$INSTDIR`` should be the path to ``TX2000``; ``$RatioBidding`` is suggested to be within 0 and 0.3; ``Timestep`` is an integer within 1 and 8760 (hours), where 1-24 hrs are tested the most;  

``$MethodSolve`` select how the problem is sloved: PSQP is solving by decomposition where Master is solved using SQP, Preg is solving by decomposition where Master is solved using IPOPT Quasi-Newton with registed user-defined functions and gradient functions; other string input leads to solving sequentially using IPOPT.

The script write files ``solutionBid_$(Timestep)h_$(RatioBidding).txt`` and ``Plantbus_$(Timestep)h_$(RatioBidding).CSV`` in the ``$OUTDIR``.