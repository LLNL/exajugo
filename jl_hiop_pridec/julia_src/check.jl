import Pkg
try
    using Glob
catch e
    if e isa ArgumentError
        println("Glob not found. Installing Glob...")
        Pkg.add("Glob")
        using Glob
    else
        rethrow(e)
    end
end

try
    using TimeZones
catch e
    if e isa ArgumentError
        println("TimeZones not found. Installing TimeZones...")
        Pkg.add("TimeZones")
        using TimeZones
    else
        rethrow(e)
    end
end

using Dates
using Glob
using TimeZones

parent_folder = "output"
exclude_dir = "scripts"
mask = "iterations*.csv"

# Check for name filter in command-line arguments
name_filter = length(ARGS) > 0 ? ARGS[1] : ""

# Check if parent folder exists
if !isdir(parent_folder)
    println("      Parent folder '$parent_folder' does not exist.")
    exit()
end

# List all subdirectories in output, excluding 'scripts'
dirs = filter(f -> isdir(f) && basename(f) != exclude_dir, joinpath.(parent_folder, readdir(parent_folder)))

# If a name is provided, filter dirs by that name
if !isempty(name_filter)
    dirs = filter(f -> occursin(name_filter, basename(f)), dirs)
end

println("")

if isempty(dirs)
    println("      No directories found in $parent_folder (excluding $exclude_dir)")
else
    # Find the most recent directory by modification time
    most_recent_dir = dirs[argmax(stat.(dirs) .|> x -> x.mtime)]
    println("      Most recent directory: ", most_recent_dir)

    # Print most recent case run
    println("      --> Most recent case run: $(basename(most_recent_dir))")

    # --- READ SLURM_JOB_ID FROM environment_vars.txt ---
    env_file = joinpath(most_recent_dir, "environment_vars.txt")
    slurm_job_id = ""

    if isfile(env_file)
        for line in eachline(env_file)
            if startswith(line, "SLURM_JOB_ID=")
                global slurm_job_id = split(line, "=", limit=2)[2]
                break
            end
        end
    end

    if !isempty(slurm_job_id)
        println("      SLURM_JOB_ID: $slurm_job_id")
    end
    # ---------------------------------------------------

    # --- NEW: Look for iterations/ subdirectory ---
    iterations_dir = joinpath(most_recent_dir, "iterations")

    if !isdir(iterations_dir)
        println("      No 'iterations' directory found in $most_recent_dir")
    else
        # Find all subdirectories with numeric names inside 'iterations'
        subdirs = filter(f -> isdir(f) && occursin(r"^\d+$", basename(f)), joinpath.(iterations_dir, readdir(iterations_dir)))
        
        # Sort numerically
        function get_number(dir)
            try
                parse(Int, basename(dir))
            catch
                -1
            end
        end
        subdirs = sort(subdirs, by=get_number)

        if isempty(subdirs)
            println("      No numeric subdirectories found in $iterations_dir")
        else
            # Get local timezone
            local_tz = localzone()

            # Collect data for table
            rows = []
            for subdir in subdirs
                num = basename(subdir)

                # Always get directory creation time (or ctime) in local time
                dir_stat = stat(subdir)
                dir_created = hasproperty(dir_stat, :birthtime) ? dir_stat.birthtime : dir_stat.ctime
                dt_create_utc = unix2datetime(dir_created)
                dt_create_local = astimezone(ZonedDateTime(dt_create_utc, tz"UTC"), local_tz)
                date_created = Dates.format(dt_create_local, "yyyy-mm-dd HH:MM:SS")

                files = glob(mask, subdir)
                if isempty(files)
                    count = 0
                    last_update = ""
                    objective = ""
                else
                    # Use most recently modified file
                    file = files[argmax(stat.(files) .|> x -> x.mtime)]
                    count = max(countlines(file) - 1, 0)

                    # Get file modification time in local time
                    file_stat = stat(file)
                    dt_update_utc = unix2datetime(file_stat.mtime)
                    dt_update_local = astimezone(ZonedDateTime(dt_update_utc, tz"UTC"), local_tz)
                    last_update = Dates.format(dt_update_local, "yyyy-mm-dd HH:MM:SS")

                    # Read objective value if more than one row
                    if count > 0
                        open(file, "r") do io
                            readline(io) # skip header
                            last_line = ""
                            for line in eachline(io)
                                last_line = line
                            end
                            if !isempty(last_line)
                                fields = split(last_line, ',')
                                if length(fields) >= 2
                                    objective = strip(fields[2])
                                else
                                    objective = ""
                                end
                            else
                                objective = ""
                            end
                        end
                    else
                        objective = ""
                    end
                end
                push!(rows, (string(num), string(count), date_created, last_update, objective))
            end

            # Calculate column widths
            headers = ["number", "# of iterations", "date created", "last update", "last objective"]
            cols = [getindex.(rows, i) for i in 1:5]
            col_widths = [maximum(length.(col)) for col in cols]
            for (i, h) in enumerate(headers)
                col_widths[i] = max(col_widths[i], length(h))
            end

            # Print header
            print("      ")
            for (h, w) in zip(headers, col_widths)
                print(rpad(h, w), "  ")
            end
            println()
            # Print separator
            print("      ")
            for w in col_widths
                print("-"^w, "  ")
            end
            println()

            # Print rows
            for row in rows
                print("      ")
                for (val, w) in zip(row, col_widths)
                    print(rpad(val, w), "  ")
                end
                println()
            end
        end
    end
end

println("")