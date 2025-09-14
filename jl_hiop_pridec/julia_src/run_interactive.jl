# save this as generate_sbatch.jl

# Check for correct number of arguments
if length(ARGS) < 1
    println("Usage: julia run_interactive.jl  case")
    exit(1)
end


function check_env_vars()
    required_vars = [
        "PATH_TO_EXAJUGO",
        "PATH_TO_INSTANCES",
        "PATH_TO_HSLLIB",
        "HIOP_INSTALL_DIR"
    ]
    missing = filter(var -> !haskey(ENV, var), required_vars)
    if isempty(missing)
        println(" --- All required environment variables are set! ---")
        return true
    else
        println(" *** Missing environment variables: ", join(missing, ", ")," ***")
        return false
    end
end

if !check_env_vars()
    println(" *** Execution aborted! ***")
    println("")
    exit()

end

instance_name = ARGS[1]
cont_file = "case"

if length(ARGS) == 2
  cont_file = ARGS[2]
end

template_dir = "./sbatch_templates"
template_file = joinpath(template_dir, "interactive.sh")

if !isfile(template_file)

    if !isfile(template_file)
        println(" *** No template file $template_file for interactive mode","! ***")
        exit(1)
    end

end

println(" --- Template file $template_file"," found! ---")

#--- read hiop submission script, but skip the first line
lines = readlines(joinpath(template_dir, "hiop.sh"))
hiop_sh = join(lines[2:end], "\n")

ENV["CONTINGENCY_FILE"]=cont_file

include("hiop.jl")

ncont = get_number_of_contingencies(instance_name)
ntasks = string(ncont+1)
ntasks = ENV["SLURM_NTASKS"]

if parse(Int, ntasks)!= parse(Int, ENV["SLURM_NTASKS"])
   println("")
   println(" # of allocated processors: ", ENV["SLURM_NTASKS"])
   println(" # of processors needed: >= ", ntasks)
   println(" *** Execution aborted! ***\n")
   exit()
end

max_iter = haskey(ENV, "MAX_ITER") ? parse(Int, ENV["MAX_ITER"]) : typemax(Int)

template = read(template_file, String)

# Replace placeholders
output = replace(template, "\$HIOP_SH" => hiop_sh)

# Replace placeholders
output = replace(output, "\$CASE" => instance_name, "\$NTASKS" => ntasks, 
            "\$CONT_FILE" => cont_file, "\$MAX_ITER" => max_iter)

output = replace(output, "srun" => "srun -n$ntasks")

insert_line = "echo \" Output directory: \$OUTPUT_DIR\"\n"
search_line = "mkdir -p \"\${OUTPUT_DIR}\"\n"

output = replace(output, search_line => search_line * insert_line)

output_dir = joinpath("output","scripts")

mkpath(output_dir)


# Write to output file
output_file = joinpath(output_dir, "run_$(instance_name).sh")

insert_line = "export BATCH_FILE=$output_file\n"
search_line = "cp \"\${BATCH_FILE}\" \"\${OUTPUT_DIR}\""

output = replace(output, search_line => insert_line * search_line)

open(output_file, "w") do io
    write(io, output)
end


println(" Script generated:")
println("")
println(output)
println("")

function yes_pressed()
    println(" Execute? ENTER=yes, any key=no")
    run(`stty raw -echo`)
    c = read(stdin, Char)
    run(`stty -raw echo`)
    println()
    return c == '\r' || c == '\n'
end

ENV["BATCH_FILE"]=output_file

chmod(output_file, 0o755)

println("Before executing, allocate an interactive debug node using the following command:")
println("   salloc -N1 -n$ntasks -ppdebug")

if yes_pressed()

  print(" --- Running command: ./$output_file  ---\n")
  run(`./$output_file`)

else

   println("\n --- Generated batch file: $output_file ---\n")
   println(" --- To run on interactive node allocation: ./$output_file ---\n\n")

end



