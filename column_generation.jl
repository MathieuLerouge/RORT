################
## LIBRAIRIES ##
################


using JuMP
using GLPK



#############
## INCLUDE ##
#############


# Include file about routes implementation
include("routes_implementation.jl")



################
## SUBPROBLEM ##
################


function solve_subproblem(dual_values,
    n, r, g, Q, es, ls, I, F,
    neighbours_in, neighbours_out, distances, times, s, q, C)

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

    # Unique route
    @constraint(model, cst_0, sum(x[1,j] for j in neighbours_out[1]) == 1)

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
        - sum(dual_values[i]*x[i,j] for i in I for j in 2:n+2
            if (i,j) in keys(distances)))


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


function solve_master(routes, n, I, distances)

    ## GLOBAL VARIABLES ##

    # Get the Number of vehicules
    K = length(routes)


    ## MODEL ##

    # Define the model
    model = Model(with_optimizer(GLPK.Optimizer))


    ## DECISION VARIABLES ##

    # Define the decision variables
    @variable(model, delta[i in 1:K], Bin)


    ## CONSTRAINTS ##

    # Delivery constraint
    # ** WIP **
    @constraint(model, cst[i in I],
        sum( sum(routes[k][i,j]
            for j in 2:n+2 if (i,j) in keys(distances))*delta[k]
        for k in 1:K) == 1)
    # @constraint(model, cst[i in I],
    #     sum( sum(routes[k][i,j]
    #         for j in 2:n+2 if (i,j) in keys(distances))*delta[k]
    #     for k in 1:K) >= 1)


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
    total_distance, nb_routes, delta_star

end


function solve_relaxed_master(routes, n, I, distances)

    ## GLOBAL VARIABLES ##

    # Get the Number of vehicules
    K = length(routes)


    ## MODEL ##

    # Define the model
    model = Model(with_optimizer(GLPK.Optimizer))


    ## DECISION VARIABLES ##

    # Define the decision variables
    @variable(model, delta[i in 1:K] >= 0)


    ## CONSTRAINTS ##

    # Deefine delivery constraint
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
function fill_graph(n_tot, I, neighbours_in, neighbours_out, distances, times)
    infinite = 10000
    added_arcs_counter = 0
    for i in I
        neighbours_in[i] = union(neighbours_in[i], 1)
        neighbours_out[i] = union(neighbours_out[i], n_tot)
        if !((1, i) in keys(distances))
            distances[(1, i)] = infinite
            times[(1, i)] = infinite
            added_arcs_counter += 1
        end
        if !((i, n_tot) in keys(distances))
            distances[(i, n_tot)] = infinite
            times[(i, n_tot)] = infinite
            added_arcs_counter += 1
        end
    end
    added_arcs_counter
end

# Function which creates initial routes of the form
# {(0, i), (i, n+1)} for i in I
function initialize_routes(n_tot, I)
    routes = []
    for i in I
        route_i = zeros(Int64, n_tot, n_tot)
        route_i[1, i] = 1
        route_i[i, n_tot] = 1
        push!(routes, route_i)
    end
    routes
end

# Function to prepare column generation
function prepare_column_generation(n_tot, I,
    neighbours_in, neighbours_out, distances, times)

    # Add arcs for initialization if needed
    added_arcs_counter = fill_graph(n_tot, I,
        neighbours_in, neighbours_out, distances, times)

    # Initialize routes for the first iteration if the master problem
    routes = initialize_routes(n_tot, I)

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
            total_distance, nb_routes, routes_selection =
                solve_master(routes, n, I, distances)
            println("Done")
            println("Number of vehicules: ", nb_routes)
            println("Total distance: ", total_distance)
            for k in 1:length(routes) if Int(routes_selection[k])==1
                route = routes[k]
                println("Selected route ", k, ":")
                display_route_as_path(route, n+2)
            end end
            println("2. Running model of the relaxed master problem")
            dual_values = solve_relaxed_master(routes, n, I, distances)
            println("Done")
            println("3. Running subproblem")
            new_route, reduced_cost = solve_subproblem(dual_values,
                n, r, g, Q, es, ls, I, F,
                neighbours_in, neighbours_out, distances, times, s, q, C)
            println("Done")
            println("Reduced cost: ", reduced_cost, "\n")
            push!(routes, new_route)
            new_routes_counter += 1
        end
    else
        while reduced_cost < -epsilon
            total_distance, nb_routes, routes_selection =
                solve_master(routes, n, I, distances)
            dual_values = solve_relaxed_master(routes, n, I, distances)
            new_route, reduced_cost = solve_subproblem(dual_values,
                n, r, g, Q, es, ls, I, F,
                neighbours_in, neighbours_out, distances, times, s, q, C)
            push!(routes, new_route)
            new_routes_counter += 1
        end
    end

    # Output
    total_distance, nb_routes, new_routes_counter

end
