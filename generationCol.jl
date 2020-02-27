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
## DATA ##
##########

# Set file_name
file_name = "E_data.txt"
#file_name = "E_data_1.txt"

# Extract data
n_tot, m, r, g, Q, es, ls, I, F,
    neighbours_in, neighbours_out, distances, times =
    extract_data(string(data_folder_path, file_name))
n = n_tot - 2

# Complete data
# Assumption: no time serive
s = [0 for i in 1:n+2]
# Assumption: no demand to deliver
q = [0 for i in 1:n+2]
# Assumption: infinite capacity of vehicules
C = 100000



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
    println("Start to run model of the sub-problem")
    optimize!(model)
    println("End of the run.")


    ## RESULTS ##

    # Get solutions and values of the objective
    x_star = value.(x)
    reduced_cost = objective_value(model)

    # Display
    println("Reduced cost: ", reduced_cost)
    println("")

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
    println("Start to run model of the master problem")
    optimize!(model)
    println("End of the run.")


    ## RESULTS ##

    # Get solutions and values of the objective
    total_distance = objective_value(model)
    delta_star = value.(delta)
    nb_routes = Int(sum(delta_star))

    # Display
    println("Current number of routes: ", nb_routes)
    println("Current total distance: ", total_distance)
    println("")

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
    println("Start to run model of the relaxed master problem")
    optimize!(model)
    println("End of the run. \n")


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
    for i in I
        neighbours_in[i] = union(neighbours_in[i], 1)
        neighbours_out[i] = union(neighbours_out[i], n+2)
        if((1,i) in keys(distances))
            continue
        else
            distances[(1, i)] = infinite
            times[(1, i)] = infinite
        end
        if((i, n+2) in keys(distances))
            continue
        else
            distances[(i, n+2)] = infinite
            times[(i, n+2)] = infinite
        end
    end
end


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



#######################
## COLUMN GENERATION ##
#######################

# Add arcs for initialization if needed
fill_graph(n_tot, m, r, g, Q, es, ls, I, F, neighbours_in, neighbours_out,
    distances, times)

# Initialize routes for the first iteration if the master problem
routes = initialize_routes()

# Initialize results with any values
global_total_distance = -1
global_nb_routes = -1

# Initialize reduced cost to a negative value (e.g. -1)
global_reduced_cost = -1

# Loop while reduced is negative
while global_reduced_cost < 0
    total_distance, nb_routes = solve_master(routes)
    dual_values = solve_relaxed_master(routes)
    new_route, reduced_cost = solve_subproblem(dual_values)
    push!(routes, new_route)
    global global_total_distance = total_distance
    global global_nb_routes = nb_routes
    global global_reduced_cost = reduced_cost
end

# Display results
println("Number of routes: ", global_nb_routes)
println("Total distance: ", global_total_distance)
