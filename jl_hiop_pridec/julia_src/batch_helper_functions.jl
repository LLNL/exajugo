

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



function yes_pressed()
    println(" Execute? ENTER=yes, any key=no")
    run(`stty raw -echo`)
    c = read(stdin, Char)
    run(`stty -raw echo`)
    println()
    return c == '\r' || c == '\n'
end

