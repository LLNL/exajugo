# save this as generate_sbatch.jl

# Check for correct number of arguments
if length(ARGS) < 2
    println("Usage: julia generate_sbatch.jl case batch_time")
    exit(1)
end

instance_name = ARGS[1]
batch_time = ARGS[2]
cont_file = "case"

if length(ARGS) == 3
  cont_file = ARGS[3]
end


# Read the template file
template_file = "template.sbatch"
if !isfile(template_file)
    println("Template file $template_file not found.")
    exit(1)
end

ENV["CONTINGENCY_FILE"]=cont_file

include("hiop.jl")

ncont = get_number_of_contingencies(instance_name)
ntasks = string(ncont+1)

max_iter = haskey(ENV, "MAX_ITER") ? parse(Int, ENV["MAX_ITER"]) : typemax(Int)


template = read(template_file, String)

# Replace placeholders
output = replace(template, "\$CASE" => instance_name, "\$NTASKS" => ntasks, 
            "\$TIME" => batch_time, "\$CONT_FILE" => cont_file, "\$MAX_ITER" => max_iter)

# Write to output file
output_file = "sub_$(instance_name).sbatch"
open(output_file, "w") do io
    write(io, output)
end

println("\n --- Generated batch file: $output_file ---\n")
println(" --- RUN: sbatch $output_file ---\n\n")

