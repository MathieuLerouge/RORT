#############
## INCLUDE ##
#############


# Set folder_path
#folder_path = "/Users/camillegrange/Documents/RORT/"
folder_path = "/Users/lerougemathieu/Documents/Courses/MPRO/RORT/Project/"

# Include file about data extraction
data_folder = "data/"
data_folder_path = string(folder_path, data_folder)
include(string(data_folder_path, "data_extraction.jl"))

# Include files implementing the different solving approaches
include("milp.jl")
include("column_generation.jl")
include("heuristic.jl")



######################
## GLOBAL VARIABLES ##
######################

## VARIABLES TO SET ##

# Set boolean for measuring running time
#measuring_time = true
measuring_time = false

# Set boolean for simplifying the graph of connections
# by removing unfeasible arcs
simplifying_graph = true
#simplifying_graph = false

# Set file_name
#file_name = "E_data.txt"
file_name = "E_data_1.txt"
#file_name = "E_data_2.txt"
#file_name = "E_data_3.txt"
#file_name = "instance_8.txt"
#file_name = "instance_10.txt"
#file_name = "instance_12.txt"
#file_name = "instance_16.txt"
#file_name = "instance_42.txt"

# Set algorithm
# - MILP: 1
# - Column generation: 2
# - Heuristic: 3
algorithm_type = 1

# Set number of duplicates for MILP
nb_duplicates = 2

# Set factor of compromise between shortest path and visit of new customers
# for heuristic
penalization_factor = 3

## OTHER VARIABLES ##

# Define boolean for displaying details of the entire solving process
# if and only if not measuring running time
display_process = !measuring_time

# Initialize time variable
total_time = 0



#####################
## DATA EXTRACTION ##
#####################

# Display
println("\n",
    "---------------\n",
    "DATA EXTRACTION\n",
    "---------------\n")

# Compile code if measuring time
# Basically, run functions once so that they are all compiled
if measuring_time
    _, _, _, _, _, _, _, _, _, _, _, _, _ =
        extract_data(string(data_folder_path, file_name), nb_duplicates)
end

# Extract data
data_extraction_time = @elapsed n_tot, m, r, g, Q, es, ls, I, F,
    neighbours_in, neighbours_out, distances, times =
        extract_data(string(data_folder_path, file_name), nb_duplicates)
n = n_tot - 2
total_time += data_extraction_time

# Complete data
# Assumption: no time service
s = [0 for i in 1:n_tot]
# Assumption: no demand to deliver
q = [0 for i in 1:n_tot]
# Assumption: infinite capacity of vehicules
C = 100000


# Compile code if measuring time
# Basically, run functions once so that they are all compiled
if measuring_time
    _ = clear_graph(n, m, r, g, Q, es, ls, I, F,
        neighbours_in, neighbours_out, distances, times, s, q, C)
end

# Delete unfeasible arcs
nb_deleted_arcs = 0
graph_simplification_time = 0
if simplifying_graph
    graph_simplification_time = @elapsed nb_deleted_arcs =
        clear_graph(n, m, r, g, Q, es, ls, I, F,
            neighbours_in, neighbours_out, distances, times, s, q, C)
    total_time += graph_simplification_time
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
    if measuring_time
        println("Elapsed time to extract data: ", data_extraction_time, "s")
    end
    if simplifying_graph
        println("Number of deleted unfeasible arcs: ", nb_deleted_arcs)
        if measuring_time
            println("Elapsed time to simplify graph: ",
                graph_simplification_time, "s")
        end
    end
end

println("\n",
    "----------------------\n",
    "END OF DATA EXTRACTION\n",
    "----------------------\n")



##########
## MILP ##
##########


# Initialize time variable
milp_solving_time = 0

if algorithm_type == 1

# Display
println("\n",
    "----\n",
    "MILP\n",
    "----\n")

# Run MILP
F_prime = F
total_distance, nb_vehicules, milp_solving_time = run_MILP(n, r, g, Q, es, ls,
    I, F_prime, neighbours_in, neighbours_out, distances, times, s, q, C)
total_time += milp_solving_time

# Display results
println("- OPTIMAL VALUES")
println("Number of vehicules: ", nb_vehicules)
println("Total distance: ", total_distance)

# Display performances
if measuring_time
    println("")
    println("- PERFORMANCES")
    println("Elapsed time to solve: ", milp_solving_time, "s")
end

println("\n",
    "-----------\n",
    "END OF MILP\n",
    "-----------\n")

end



#######################
## COLUMN GENERATION ##
#######################

# Initialize time variables
initialization_time = 0
column_generation_time = 0

if algorithm_type == 2

# Display
println("\n",
    "-----------------\n",
    "COLUMN GENERATION\n",
    "-----------------\n")

# Compile code if measuring time
# Basically, run functions once so that they are all compiled
if measuring_time
    routes, _ = prepare_column_generation(n_tot, I,
        neighbours_in, neighbours_out, distances, times)
    _, _, _ = run_column_generation(routes, false)
end

# Prepare column generation
# In particular, initialize a first set of routes
initialization_time = @elapsed routes, added_arcs_counter =
    prepare_column_generation(n_tot, I,
        neighbours_in, neighbours_out, distances, times)
total_time += initialization_time
println("- GRAPH ALTERATIONS")
println("Number of infinite arcs added to the graph: ", added_arcs_counter)
println("")

# Display initial routes
route_index = 1
if display_process
    println("- INITIALIZATION")
    println("Number of initial routes: ", length(routes))
    for route in routes
        println("Initial route ", route_index, ":")
        display_route_as_path(route, n_tot)
        global route_index += 1
    end
    println("")
end

# Run column generation
column_generation_time = @elapsed total_distance, nb_routes,
    new_routes_counter = run_column_generation(routes, display_process)
total_time += column_generation_time

# Display results
println("- OPTIMAL VALUES")
println("Number of vehicules: ", nb_routes)
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

end



###############
## HEURISTIC ##
###############

# Initialize time and performances variables
heuristic_time = 0
nb_milps = 0

if algorithm_type==3

# Display
println("\n",
    "---------\n",
    "HEURISTIC\n",
    "---------\n")

# Run heuristic
heuristic_time = @elapsed total_distance, nb_routes, routes,
    nb_milps = run_heuristic(penalization_factor,
        n, m, r, g, Q, es, ls, I, F,
        neighbours_in, neighbours_out, distances, times,
        display_process)
total_time += heuristic_time

# Display results
println("- OPTIMAL VALUES")
println("Number of vehicules: ", nb_routes)
println("Total distance: ", total_distance)
println("")

# Display performances
println("- PERFORMANCES")
println("Number of MILPS solved: ", nb_milps)
if measuring_time
    println("Elapsed time to run heuristic: ", heuristic_time, "s")
end

# Display
println("\n",
    "----------------\n",
    "END OF HEURISTIC\n",
    "----------------\n")

end



##########################
## PERFORMANCES SUMMARY ##
##########################


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
    if algorithm_type == 2
        println("Number of infinite arcs added to the graph: ",
            added_arcs_counter)
        println("Number of routes added / of subproblems solved: ",
            new_routes_counter)
    elseif algorithm_type == 3
        println("Elapsed time to run heuristic: ",
            heuristic_time, "s")
    end
    println("")

    println("- TIMES")
    println("Elapsed time to extract data: ", data_extraction_time, "s")
    if simplifying_graph
        println("Elapsed time to simplify graph: ", graph_simplification_time, "s")
    end
    if algorithm_type == 1
        println("Elapsed time to solve: ", milp_solving_time, "s")
    elseif algorithm_type == 2
        println("Elapsed time to prepare column generation: ",
            initialization_time, "s")
        println("Elapsed time to run column generation: ",
            column_generation_time, "s")
    elseif algorithm_type == 3
        println("Elapsed time to run heuristic: ", heuristic_time, "s")
    end
    println("Elapsed total time: ", total_time)

    println("\n",
        "--------------\n",
        "END OF SUMMARY\n",
        "--------------\n")

end
