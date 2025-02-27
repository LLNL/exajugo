#!/bin/sh
#SBATCH --job-name=testmj1mj2
#SBATCH --nodes=1
#SBATCH --tasks-per-node=18
#SBATCH --partition=pdebug
#SBATCH --time=00:30:00

### Network 1, scenario 1

printf "Running MyJulia1 ..."
srun julia -e 'include("../../MyJulia1.jl"); MyJulia1("../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_01R-10/scenario_1/case.con", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_01R-10/case.inl", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_01R-10/scenario_1/case.raw", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_01R-10/case.rop", 600, 1, "Network_01R-10")'
printf "\n\n\n"

printf "Running MyJulia2 ..."
srun julia -e 'include("../../MyJulia2.jl"); MyJulia2("../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_01R-10/scenario_1/case.con", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_01R-10/case.inl", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_01R-10/scenario_1/case.raw", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_01R-10/case.rop", 900, 1, "Network_01R-10")'
printf "\n\n\n"

printf "Running GoCompetition evaluation script ..."
python ../../../PyEvaluation/test.py "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_01R-10/scenario_1/case.raw" "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_01R-10/case.rop" "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_01R-10/scenario_1/case.con" "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_01R-10/case.inl" solution1.txt solution2.txt summary.csv detail.csv

### Network 10, scenario 1

#printf "Running MyJulia1 ..."
#srun julia -e 'include("../../MyJulia1.jl"); MyJulia1("../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_1/case.con", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_1/case.inl", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_1/case.raw", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/case.rop", 2400, 1, "Network_10R-10")'
#printf "\n\n\n"
#
#printf "Running MyJulia2 ..."
#srun julia -e 'include("../../MyJulia2.jl"); MyJulia2("../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_1/case.con", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_1/case.inl", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_1/case.raw", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/case.rop", 28000, 1, "Network_10R-10")'
#printf "\n\n\n"
#
#printf "Running GoCompetition evaluation script ..."
#python ../../../PyEvaluation/test.py "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_1/case.raw" "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/case.rop" "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_1/case.con" "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_1/case.inl" solution1.txt solution2.txt summary.csv detail.csv

### Network 10, scenario 8

#printf "Running MyJulia1 ..."
#srun julia -e 'include("../../MyJulia1.jl"); MyJulia1("../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_8/case.con", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_8/case.inl", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_8/case.raw", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/case.rop", 2400, 1, "Network_10R-10")'
#printf "\n\n\n"
#
#printf "Running MyJulia2 ..."
#srun julia -e 'include("../../MyJulia2.jl"); MyJulia2("../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_8/case.con", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_8/case.inl", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_8/case.raw", "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/case.rop", 28000, 1, "Network_10R-10")'
#printf "\n\n\n"
#
#printf "Running GoCompetition evaluation script ..."
#python ../../../PyEvaluation/test.py "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_8/case.raw" "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/case.rop" "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_8/case.con" "../../../goinstances/challenge1/Original_Dataset_1-4/Original_Dataset_Real-Time_Edition_1/Network_10R-10/scenario_8/case.inl" solution1.txt solution2.txt summary.csv detail.csv
