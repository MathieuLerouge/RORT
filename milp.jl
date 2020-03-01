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
#file_name = "E_data_2.txt"
file_name = "E_data_3.txt"


## DATA EXTRACTION ##

# Display
println("\n",
    "---------------\n",
    "DATA EXTRACTION\n",
    "---------------\n")

# Compile code if measuring time
# Basically, run functions once so that they are all compiled
if measuring_time
    n_tot, m, r, g, Q, es, ls, I, F_prime, neighbours_in, neighbours_out,
        distances, times = extract_data(string(data_folder_path, file_name))
    n = n_tot - 2
    _ = clear_graph(n, m, r, g, Q, es, ls, I, F_prime,
        neighbours_in, neighbours_out, distances, times, s, q, C)
end

# Extract data
data_extraction_time = @elapsed n_tot, m, r, g, Q, es, ls, I, F_prime,
    neighbours_in, neighbours_out, distances, times =
        extract_data(string(data_folder_path, file_name))
n = n_tot - 2

# Complete data
# Assumption: no time service
s = [0 for i in 1:n+2]
# Assumption: no demand to deliver
q = [0 for i in 1:n+2]
# Assumption: infinite capacity of vehicules
C = 100000

# Delete unfeasible arcs
nb_deleted_arcs = 0
graph_simplification_time = 0
if simplifying_graph
    graph_simplification_time = @elapsed nb_deleted_arcs =
        clear_graph(n, m, r, g, Q, es, ls, I, F_prime,
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


## MILP ##

# Display
println("\n",
    "----\n",
    "MILP\n",
    "----\n")

# Define the model
model = Model(with_optimizer(GLPK.Optimizer))

# Define the decision variables
@variable(model, x[i in 1:n+1,j in 2:n+2], Bin)
@variable(model, tau[i in 1:n+2] >= 0)
# Assumption: no demand to deliver
#@variable(model, u[i in 1:n+2] >= 0)
@variable(model, y[i in 1:n+2] >= 0)

# Define the constraints

# - Delivery constraint (3)
@constraint(model, ct_2[i in I],
    sum(x[i,j] for j in neighbours_out[i] if j!=1) == 1)

# - Recharging station constraint (3)
@constraint(model, ct_3[i in F_prime],
    sum(x[i,j] for j in neighbours_out[i] if j!=1) <= 1)

# - Flow conservation (4)
@constraint(model, ct_4[j in 2:n+1],
    sum(x[j,i] for i in neighbours_out[j] if i != 1) -
    sum(x[i,j] for i in neighbours_in[j] if i != n+2) == 0)

# - Time feasability when leaving customers and depot (5)
@constraint(model, ct_5[i in union(I, 1), j in 2:n+2 ; (i,j) in keys(distances)],
    tau[i]+(times[(i,j)] + s[i])*x[i,j] - ls[1]*(1-x[i,j]) <= tau[j])

# - Time feasability when leaving recharging stations (6)
@constraint(model, ct_6[i in F_prime, j in 2:n+2;
    (i,j) in keys(distances)],
    tau[i]+times[(i,j)]*x[i,j] + g*(Q-y[i]) - (ls[1] + g*Q)*(1-x[i,j]) <= tau[j])

# - Time windows (7)
@constraint(model, ct_7[j in 1:n+2], es[j] <= tau[j])
@constraint(model, ct_7bis[j in 1:n+2], tau[j] <= ls[j])

# - Demand feasability (8) and (9)
# Assumption: no demand to deliver
#@constraint(model, ct_8[i in 1:n+1, j in 2:n+2 ; (i,j) in keys(distances)],
#    u[j] <= u[i] - q[i]*x[i,j] + C*(1-x[i,j]))
#@constraint(model, ct_9, u[1] <= C)

# - Battery feasability (10) and (11)
@constraint(model, ct_10[i in I, j in 2:n+2;
    (i,j) in keys(distances)],
    y[j] <= y[i] - r*distances[(i,j)]*x[i,j] + Q*(1-x[i,j]))
@constraint(model, ct_11[i in union(F_prime,1), j in 2:n+2;
    (i,j) in keys(distances) && i != n+2],
    y[j] <= Q - r*distances[(i,j)]*x[i,j])

# Define the objective function
@objective(model, Min,
    sum( distances[(i,j)]*x[i,j]
        for (i,j) in keys(distances) if i != n+2 && j!= 1))

# Run the solving process
solving_time = @elapsed optimize!(model)

# Get solutions
xs_star = value.(x)
nb_vehicules = Int(sum(xs_star[1,:]))

# Display results
println("- OPTIMAL VALUES")
println("Number of vehicules: ", nb_vehicules)
println("Total distance: ", objective_value(model))

# Display performances
if measuring_time
    println("")
    println("- PERFORMANCES")
    println("Elapsed time to solve: ", solving_time, "s")
end

println("\n",
    "-----------\n",
    "END OF MILP\n",
    "-----------\n")


## PERFORMANCES ##

# Display
if measuring_time

    println("\n",
        "-------\n",
        "SUMMARY\n",
        "-------\n")

    if simplifying_graph
        println("- NUMBERS")
        println("Number of deleted unfeasible arcs: ", nb_deleted_arcs)
    end

    println("- TIMES")
    total_time = data_extraction_time + graph_simplification_time +
        solving_time
    println("Elapsed time to extract data: ", data_extraction_time, "s")
    if simplifying_graph
        println("Elapsed time to simplify graph: ",
            graph_simplification_time, "s")
    end
    println("Elapsed time to solve: ", solving_time, "s")
    println("Elapsed total time: ", total_time)

    println("\n",
        "--------------\n",
        "END OF SUMMARY\n",
        "--------------\n")

end
