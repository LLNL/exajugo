# read command line parameters

INSTDIR = ARGS[1]

# load modules

println("\nLoading modules ...")
clock = @elapsed begin
	# internal modules
	push!(LOAD_PATH, joinpath(dirname(PROGRAM_FILE), "../modules"))
	using InstanceReader
end
println("All modules loaded in ", round(clock, digits=2), " secs.")
flush(stdout)

# run test over all instances within INSTDIR

numsuccess = 0
numfailed = 0
faildirs = String[]
println("\nTesting parser for all instances within ", INSTDIR, ".")
for (root, dirs, files) in walkdir(INSTDIR)
	for file in files
		if file == "case.raw"
			print("Testing for instance at ", root, " ... ")
			try
				ParseInstance(root);
				global numsuccess += 1;
				println("success!")
			catch
				push!(faildirs, root);
				global numfailed += 1;
			end
		end
	end
end

# print results

println("\nSuccess: ", numsuccess)
println("Failed:  ", numfailed)
if numfailed > 0
	println("Failed instances:")
	for dir in faildirs
		println("\t", dir)
	end
end
