# read command line parameters

MYEXE1 = ARGS[1]
INSTDIR = ARGS[2]
OUTDIR = ARGS[3]

# load modules

println("\nLoading modules ...")
clock = @elapsed begin
	# internal modules
	push!(LOAD_PATH, joinpath(dirname(PROGRAM_FILE), "../modules"))
	import InstanceReader: instancefilenames
end
println("All modules loaded in ", round(clock, digits=2), " secs.")
flush(stdout)

# function to create execution script

function preparerun(instancedir::String, instanceid::Int)::Nothing
	exedir = OUTDIR*"/testrun_"*string(instanceid)
	mkpath(exedir)
	f = open(exedir*"/submit.sh", "w")
	println(f, "#!/bin/sh\n",
		"#SBATCH --account=gridopt\n",
		"#SBATCH --nodes=1\n",
		"#SBATCH --tasks-per-node=18\n",
		"#SBATCH --partition=pdebug\n",
		"#SBATCH --time=00:15:00\n",
		"#SBATCH --job-name=testME1_", instanceid, "\n\n",
		"mpiexec -n 18 ", MYEXE1, " case.con case.inl case.raw case.rop 600 1 SomeNetwork\n")
	close(f)
	files = instancefilenames(instancedir)
	for f in files
		symlink(f, exedir*"/case."*f[(end-2):end])
	end
end

# create execution scripts for all instances within INSTDIR

println("\nPreparing submission scripts for all instances within ", INSTDIR, ".")
id = 1
RootDirectories = String[]
for (root, dirs, files) in walkdir(INSTDIR)
	for file in files
		if file == "case.raw"
			println("Creating execution folder and submission directory for instance at ", root, ".")
			push!(RootDirectories, root)
			preparerun(root, id)
			global id += 1
		end
	end
end

# create master launcher and register of execution scripts

fmaster = open(OUTDIR*"/launch_all_tests.sh", "w")
println(fmaster, "#!/bin/bash\n")
freg = open(OUTDIR*"/test_register.csv", "w")
println(freg, "id,root")
for i = 1:(id-1)
	println(fmaster, "cd testrun_", i, " && sbatch submit.sh && cd ..")
	println(freg, i, ",", RootDirectories[i])
end
close(fmaster)
close(freg)
