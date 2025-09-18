# save this as generate_sbatch.jl

include("batch_helper_functions.jl")

# Check for correct number of arguments
if length(ARGS) < 2
    println("Usage: julia generate_sbatch.jl case batch_time")
    exit(1)
end


if !check_env_vars()
    println(" *** Execution aborted! ***")
    println("")
    exit()

end

instance_name = ARGS[1]
batch_time = ARGS[2]
cont_file = "case"

if length(ARGS) == 3
  cont_file = ARGS[3]
end

template_dir = "./sbatch_templates"
template_file = joinpath(template_dir, "default.sbatch")

if !isfile(template_file)

    # Get hostname
    hostname = gethostname()

    # Remove digits from hostname
    base_hostname = replace(hostname, r"\d+" => "")

    template_file = joinpath(template_dir, base_hostname * ".sbatch")

    if !isfile(template_file)
        println(" *** No template file $template_file","! ***")
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
ntasks = ncont+1
if haskey(ENV, "NTASKS")
   ntasks = parse(Int, ENV["NTASKS"])
else
  new_value = read_ntasks_or_nothing()
  if new_value !== nothing
     ntasks = new_value
  end
end

error_tasks(ntasks)
ntasks = string(ntasks)

max_iter = haskey(ENV, "MAX_ITER") ? parse(Int, ENV["MAX_ITER"]) : typemax(Int)

template = read(template_file, String)

# Replace placeholders
output = replace(template, "\$HIOP_SH" => hiop_sh)

# Replace placeholders
output = replace(output, "\$CASE" => instance_name, "\$NTASKS" => ntasks, 
            "\$TIME" => batch_time, "\$CONT_FILE" => cont_file, "\$MAX_ITER" => max_iter)

output_dir = joinpath("output","scripts")

mkpath(output_dir)

# Write to output file
output_file = joinpath(output_dir, "sub_$(instance_name).sbatch")
open(output_file, "w") do io
    write(io, output)
end

println(" Script generated:")
println("")
println(output)
println("")

ENV["BATCH_FILE"]=output_file

if yes_pressed()

   run(`sbatch $output_file`)
   println("\n --- Batch submitted: $output_file ---\n")


else

   println("\n --- Generated batch file: $output_file ---\n")
   println(" --- To submit the job, run: sbatch $output_file ---\n\n")

end



