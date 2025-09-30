

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


function read_ntasks_or_nothing()
    print(" Enter # of tasks or press ENTER to set # of tasks= # of contingencies+1): ")
    input = readline()
    if isempty(input)
        return nothing
    else
        try
            return parse(Int, input)
        catch e
            println("Invalid input. Please enter a valid integer or press Enter.")
            return read_ntasks_or_nothing()
        end
    end
end


function error_tasks(ntasks)

   if ntasks<1
       println(" *** # of tasks must be >= 1 ***")
       exit(1)
   end

end

function print_generated_script(output)
    prn_env = get(ENV, "PRN_SCRIPT", nothing)
    if prn_env === nothing
        # Environment variable does not exist; print output
        println(" Script generated:")
        println("")
        println(output)
        println("")
    else
        prn_env_lower = lowercase(prn_env)
        if prn_env_lower == "true" || prn_env_lower == "yes" || prn_env_lower == "1"
            println(" Script generated:")
            println("")
            println(output)
            println("")
        end
        # If value is "false", "no", or "0", do nothing
    end
end


function yes_pressed()
    exec_env = get(ENV, "EXECUTE", nothing)
    if exec_env !== nothing
        exec_env_lower = lowercase(exec_env)
        if exec_env_lower == "yes" || exec_env_lower == "1" || exec_env_lower == "true"
            return true
        else
            return false
        end
    end

    println(" Execute? ENTER=yes, any key=no")
    run(`stty raw -echo`)
    c = read(stdin, Char)
    run(`stty -raw echo`)
    println()
    return c == '\r' || c == '\n'
end


