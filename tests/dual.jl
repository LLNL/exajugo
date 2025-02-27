using JuMP, GLPK, Random, DataFrames

# Define a function to create constraints
function add_constraints(model::Model, x, y, n::Int, k::Int, constraint_status)
    for i in 1:n
        for j in 1:n
            # @constraint(model, Pconst[i,j,k], x[i] + y[j] >= i + j + k)
            @constraint(model, x[i] + y[j] >= i + j + k)
            push!(constraint_status, (i=i, j=j, k=k))
        end
    end
end

# Define the main model function
function main_model(n::Int, iterations::Int)
    model = Model(GLPK.Optimizer)
    
    # Define variables
    @variable(model, x[1:n] >= 0)
    @variable(model, y[1:n] >= 0)
    
    # Data structure to store constraint statuses
    constraint_status = []
    
    for k in 1:iterations
        # Add constraints for the current iteration
        add_constraints(model, x, y, n, k, constraint_status)
    end
    
    # Add objective function or additional constraints
    
    # Solve the model
    optimize!(model)
    
    # Initialize DataFrame columns
    df = DataFrame(i = Int[], j = Int[], k = Int[], status = Int[])  # Change status type to Int
    
    # Get status of constraints after optimization and save to DataFrame
    for (i, j, k) in constraint_status
        constraint_value = value(x[i] + y[j] - (i + j + k))
        status = isapprox(constraint_value, 0, atol=1e-6) ? 1 : 0  # Change status to 1 for active, 0 for inactive
        push!(df, (i=i, j=j, k=k, status=status))
    end
    
    return df
end

# Generate random data
Random.seed!(123)  # for reproducibility
n = 3
iterations = 2  # Just for demonstration

# Call the main model function
df = main_model(n, iterations)
println(df)
