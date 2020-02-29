################
## LIBRAIRIES ##
################

using JuMP
using GLPK



#############
## INCLUDE ##
#############

# Set folder_path
#folder_path = "/Users/camillegrange/Documents/RORT/"
folder_path = "/Users/lerougemathieu/Documents/Courses/MPRO/RORT/Project/"

# Prepare data extraction
data_folder = "data/"
data_folder_path = string(folder_path, data_folder)
include(string(data_folder_path, "data_extraction.jl"))



################
## SUBPROBLEM ##
################

function solve_subproblem(dual_values)

    ## MODEL ##

    # Define the model
    model = Model(with_optimizer(GLPK.Optimizer))


    ## DECISION VARIABLES ##

    # Define the decision variables
    @variable(model, x[i in 1:n+1,j in 2:n+2], Bin)
    @variable(model, tau[i in 1:n+2] >= 0)
    # Assumption: no demand to deliver
    #@variable(model, u[i in 1:n+2] >= 0)
    @variable(model, y[i in 1:n+2] >= 0)


    ## CONSTRAINTS ##

    # Flow conservation (4)
    @constraint(model, cst_4[j in 2:n+1],
        sum(x[j,i] for i in neighbours_out[j]) -
        sum(x[i,j] for i in neighbours_in[j]) == 0)

    # Time feasability when leaving customers and depot (5)
    @constraint(model, cst_5[i in union(I, 1), j in 2:n+2; (i,j) in keys(distances)],
        tau[i]+(times[(i,j)] + s[i])*x[i,j] - ls[1]*(1-x[i,j]) <= tau[j])

    # Time feasability when leaving recharging stations (6)
    @constraint(model, cst_6[i in F, j in 2:n+2; (i,j) in keys(distances)],
        tau[i]+times[(i,j)]*x[i,j] + g*(Q-y[i]) - (ls[1] + g*Q)*(1-x[i,j])
        <= tau[j])

    # Time windows (7)
    @constraint(model, cst_7[j in 1:n+2], es[j] <= tau[j])
    @constraint(model, cst_7bis[j in 1:n+2], tau[j] <= ls[j])

    # Demand feasability (8) and (9)
    # Assumption: no demand to deliver
    #@constraint(model, cst_8[i in 1:n+1, j in 2:n+2; (i,j) in keys(distances)],
    #    u[j] <= u[i] - q[i]*x[i,j] + C*(1-x[i,j]))
    #@constraint(model, cst_9, u[1] <= C)

    # Battery feasability (10) and (11)
    @constraint(model, cst_10[i in I, j in 2:n+2;
        (i,j) in keys(distances)],
        y[j] <= y[i] - r*distances[(i,j)]*x[i,j] + Q*(1-x[i,j]))
    @constraint(model, cst_11[i in union(F,1), j in 2:n+2;
        (i,j) in keys(distances)],
        y[j] <= Q - r*distances[(i,j)]*x[i,j])


    ## OBJECTIVE ##

    # Define the objective function
    @objective(model, Min,
        sum(distances[(i,j)]*x[i,j] for (i,j) in keys(distances))
        - sum(dual_values[i]*x[i,j] for i in I for j in 2:n+2 if (i,j) in keys(distances)))


    ## SOLVING ##

    # Run the solving process
    optimize!(model)


    ## RESULTS ##

    # Get solutions and values of the objective
    x_star = value.(x)
    reduced_cost = objective_value(model)

    # Output
    x_star, reduced_cost

end



####################
## MASTER PROBLEM ##
####################

function solve_master(routes)

    ## GLOBAL VARIABLES ##

    # Get the number of routes
    K = length(routes)


    ## MODEL ##

    # Define the model
    model = Model(with_optimizer(GLPK.Optimizer))


    ## DECISION VARIABLES ##

    # Define the decision variables
    @variable(model, delta[i in 1:K], Bin)


    ## CONSTRAINTS ##

    # Delivery constraint
    @constraint(model, cst[i in I],
        sum( sum(routes[k][i,j]
            for j in 2:n+2 if (i,j) in keys(distances))*delta[k]
        for k in 1:K) == 1)


    ## OBJECTIVE ##

    # Define the objective function
    @objective(model, Min,
        sum( sum(routes[k][i,j]*distances[(i,j)]
            for i in 1:n+1 for j in 2:n+2 if (i,j) in keys(distances))*delta[k]
        for k in 1:K))


    ## SOLVING ##

    # Run the solving process
    optimize!(model)


    ## RESULTS ##

    # Get solutions and values of the objective
    total_distance = objective_value(model)
    delta_star = value.(delta)
    nb_routes = Int(sum(delta_star))

    # Output
    total_distance, nb_routes

end


function solve_relaxed_master(routes)

    ## GLOBAL VARIABLES ##

    # Get the number of routes
    K = length(routes)


    ## MODEL ##

    # Define the model
    model = Model(with_optimizer(GLPK.Optimizer))


    ## DECISION VARIABLES ##

    # Define the decision variables
    @variable(model, delta[i in 1:K] >= 0)


    ## CONSTRAINTS ##

    # Delivery constraint
    @constraint(model, cst_bin[i in K], delta[i] <= 1)
    @constraint(model, cst[i in I],
        sum( sum(routes[k][i,j]
            for j in 2:n+2 if (i,j) in keys(distances))*delta[k]
        for k in 1:K) == 1)


    ## OBJECTIVE ##

    # Define the objective function
    @objective(model, Min,
        sum( sum(routes[k][i,j]*distances[(i,j)]
            for i in 1:n+1 for j in 2:n+2 if (i,j) in keys(distances))*delta[k]
        for k in 1:K))


    ## SOLVING ##

    # Run the solving process
    optimize!(model)


    ## RESULTS ##

    # Get values of duals variables
    dual_values = Dict{Int64,Float64}()
    for i in I
        dual_values[i] = dual(cst[i])
    end

    # Output
    dual_values

end



######################################
## INITIALIZATION OF MASTER PROBLEM ##
######################################

# Function which adds arcs (0, i) and (i, n+1)
# with infinite distances and times, if needed
function fill_graph(n_tot, m, r, g, Q, es, ls, I, F,
    neighbours_in, neighbours_out, distances, times)
    infinite = 10000
    added_arcs_counter = 0
    for i in I
        neighbours_in[i] = union(neighbours_in[i], 1)
        neighbours_out[i] = union(neighbours_out[i], n+2)
        if !((1,i) in keys(distances))
            distances[(1, i)] = infinite
            times[(1, i)] = infinite
            added_arcs_counter += 1
        end
        if !((i, n+2) in keys(distances))
            distances[(i, n+2)] = infinite
            times[(i, n+2)] = infinite
            added_arcs_counter += 1
        end
    end
    added_arcs_counter
end

# Function which creates initial routes of the form
# {(0, i), (i, n+1)} for i in I
function initialize_routes()
    routes = []
    for i in I
        route_i = zeros(Int64, n+2, n+2)
        route_i[1, i] = 1
        route_i[i, n+2] = 1
        push!(routes, route_i)
    end
    routes
end

# Function to prepare column generation
function prepare_column_generation(n_tot, m, r, g, Q, es, ls, I, F,
    neighbours_in, neighbours_out, distances, times)

    # Add arcs for initialization if needed
    added_arcs_counter = fill_graph(n_tot, m, r, g, Q, es, ls, I, F,
        neighbours_in, neighbours_out, distances, times)

    # Initialize routes for the first iteration if the master problem
    routes = initialize_routes()

    # Output
    routes, added_arcs_counter

end



#######################
## COLUMN GENERATION ##
#######################

function run_column_generation(routes, display_process)

    # Set accuracy
    epsilon = 0.000001

    # Initialize variables storing results with any initial values
    total_distance = -1
    nb_routes = -1

    # Initialize variable storing reduced cost to a negative value (e.g. -1)
    reduced_cost = -1

    # Initialize variable storing the number of subproblem solved to 0
    new_routes_counter = 0

    # Loop while reduced is negative
    if display_process
        while reduced_cost < -epsilon
            println("- ITERATION ", new_routes_counter)
            println("1. Running model of the master problem")
            total_distance, nb_routes = solve_master(routes)
            println("Done")
            println("Number of routes: ", nb_routes)
            println("Total distance: ", total_distance)
            println("2. Running model of the relaxed master problem")
            dual_values = solve_relaxed_master(routes)
            println("Done")
            println("3. Running subproblem")
            new_route, reduced_cost = solve_subproblem(dual_values)
            println("Done")
            println("Reduced cost: ", reduced_cost, "\n")
            push!(routes, new_route)
            new_routes_counter += 1
        end
    else
        while reduced_cost < -epsilon
            total_distance, nb_routes = solve_master(routes)
            dual_values = solve_relaxed_master(routes)
            new_route, reduced_cost = solve_subproblem(dual_values)
            push!(routes, new_route)
            new_routes_counter += 1
        end
    end

    # Output
    total_distance, nb_routes, new_routes_counter

end



##########
## MAIN ##
##########


## VARIABLES TO SET ##

# Set boolean for measuring running time
measuring_time = true
#measuring_time = false

# Set boolean for simplifying the graph of connections
# by removing unfeasible arcs
simplifying_graph = true
#simplifying_graph = false

# Set file_name
#file_name = "E_data.txt"
#file_name = "E_data_1.txt"
file_name = "E_data_2.txt"
#file_name = "E_data_3.txt"


## OTHER VARIABLES ##

# Define boolean for displaying details of the entire solving process
# if and only if not measuring running time
display_process = !measuring_time


## DATA EXTRACTION ##

# Display
println("\n",
    "---------------\n",
    "DATA EXTRACTION\n",
    "---------------\n")

# Compile code if measuring time
# Basically, run functions once so that they are all compiled
if measuring_time
    n_tot, m, r, g, Q, es, ls, I, F, neighbours_in, neighbours_out,
        distances, times = extract_data(string(data_folder_path, file_name))
    n = n_tot - 2
    _ = clear_graph(n, m, r, g, Q, es, ls, I, F,
        neighbours_in, neighbours_out, distances, times, s, q, C)
end

# Extract data
data_extraction_time = @elapsed n_tot, m, r, g, Q, es, ls, I, F,
    neighbours_in, neighbours_out, distances, times =
        extract_data(string(data_folder_path, file_name))
n = n_tot - 2

# Complete data
# Assumption: no time service
s = [0 for i in 1:n+2]
# Assumption: no demand to deliver
#q = [0 for i in 1:n+2]
# Assumption: infinite capacity of vehicules
#C = 100000

# Delete unfeasible arcs
nb_deleted_arcs = 0
graph_simplification_time = 0
if simplifying_graph
    graph_simplification_time = @elapsed nb_deleted_arcs =
        clear_graph(n, m, r, g, Q, es, ls, I, F,
            neighbours_in, neighbours_out, distances, times, s, q, C)
end

# Display data
if display_process
    println("- DATA")
    print("Indices of clients I: ")
    for i in I
        print(i, ",")
    end
    println("")
    print("Indices of stations F: ")
    for i in F
        print(i, ",")
    end
    if simplifying_graph
        println("\n")
    else
        println("")
    end
end

# Display performances
if simplifying_graph || measuring_time
    println("- PERFORMANCES")
    if simplifying_graph
        println("Number of deleted unfeasible arcs: ", nb_deleted_arcs)
        println("Elapsed time to simplify graph: ", graph_simplification_time, "s")
    end
    if measuring_time
        println("Elapsed time to extract data: ", data_extraction_time, "s")
    end
end

println("\n",
    "----------------------\n",
    "END OF DATA EXTRACTION\n",
    "----------------------\n")


## COLUMN GENERATION ##

# Display
println("\n",
    "-----------------\n",
    "COLUMN GENERATION\n",
    "-----------------\n")

# Compile code if measuring time
# Basically, run functions once so that they are all compiled
if measuring_time
    _, _ = prepare_column_generation(n_tot, m, r, g, Q, es, ls, I, F,
        neighbours_in, neighbours_out, deepcopy(distances), deepcopy(times))
    _, _, _ = run_column_generation(routes, false)
end

# Prepare column generation
# In particular, initialize a first set of routes
initialization_time = @elapsed routes, added_arcs_counter =
    prepare_column_generation(n_tot, m, r, g, Q, es, ls, I, F,
        neighbours_in, neighbours_out, distances, times)
println("- GRAPH ALTERATIONS")
println("Number of infinite arcs added to the graph: ", added_arcs_counter)
println("")

# Display initial routes
route_index = 1
if display_process
    println("- INITIALIZATION")
    println("Number of initial routes: ", length(routes))
    for route in routes
        println("Route ", route_index, ":")
        for i in 1:n+2
            println(route[i, :])
        end
        global route_index += 1
    end
    println("")
end

# Run column generation
column_generation_time = @elapsed total_distance, nb_routes,
    new_routes_counter = run_column_generation(routes, display_process)

# Display results
println("- OPTIMAL VALUES")
println("Number of routes: ", nb_routes)
println("Total distance: ", total_distance)
println("")

# Display performances
println("- PERFORMANCES")
println("Number of routes added / of subproblems solved: ", new_routes_counter)
if measuring_time
    println("Elapsed time to prepare column generation: ",
        initialization_time, "s")
    println("Elapsed time to run column generation: ",
        column_generation_time, "s")
end

println("\n",
    "------------------------\n",
    "END OF COLUMN GENERATION\n",
    "------------------------\n")



## PERFORMANCES ##

# Display
if measuring_time

    println("\n",
        "-------\n",
        "SUMMARY\n",
        "-------\n")

    println("- NUMBERS")
    if simplifying_graph
        println("Number of deleted unfeasible arcs: ", nb_deleted_arcs)
    end
    println("Number of infinite arcs added to the graph: ", added_arcs_counter)
    println("Number of routes added / of subproblems solved: ", new_routes_counter)
    println("")

    println("- TIMES")
    total_time = data_extraction_time + graph_simplification_time +
        initialization_time + column_generation_time
    println("Elapsed time to extract data: ", data_extraction_time, "s")
    if simplifying_graph
        println("Elapsed time to simplify graph: ", graph_simplification_time, "s")
    end
    println("Elapsed time to prepare column generation: ",
        initialization_time, "s")
    println("Elapsed time to run column generation: ",
        column_generation_time, "s")
    println("Elapsed total time: ", total_time)

    println("\n",
        "--------------\n",
        "END OF SUMMARY\n",
        "--------------\n")

end
