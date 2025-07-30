# save this as generate_sbatch.jl

# Check for correct number of arguments
if length(ARGS) < 2
    println("Usage: julia generate_sbatch.jl case batch_time")
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

check_env_vars()

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
ntasks = string(ncont+1)

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

function yes_pressed()
    println(" Execute? ENTER=yes, any key=no")
    run(`stty raw -echo`)
    c = read(stdin, Char)
    run(`stty -raw echo`)
    println()
    return c == '\r' || c == '\n'
end

if yes_pressed()

  run(`sbatch $output_file`)

else

   println("\n --- Generated batch file: $output_file ---\n")
   println(" --- To submit the job, run: sbatch $output_file ---\n\n")

end



