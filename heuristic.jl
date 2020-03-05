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



##########################
## GENERATION OF ROUTES ##
##########################


function generate_route(must_be_visited_I,
    non_visited_I, avg_distance, penalization_factor,
    n, r, g, Q, es, ls, I, F,
    neighbours_in, neighbours_out, distances, times)

    ## VARIABLES ##
    n_tot = n+2

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
    @constraint(model, cst_unique, sum(x[1,j] for j in neighbours_out[1]) == 1)

    # Compulsory visit
    @constraint(model, cst_visit[i in must_be_visited_I],
        sum(x[i,j] for j in neighbours_out[i]) == 1)

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
    if (penalization_factor <= 0)
        @objective(model, Min,
            sum(distances[(i,j)]*x[i,j] for (i,j) in keys(distances)))
    else
        @objective(model, Min,
            sum(distances[(i,j)]*x[i,j] for (i,j) in keys(distances))
            - penalization_factor*avg_distance*sum(
            x[i,j] for (i,j) in keys(distances) if i in non_visited_I))
    end


    ## SOLVING ##

    # Run the solving process
    optimize!(model)


    ## RESULTS ##

    # Get solution
    x_star = value.(x)

    # Compute elements of the route and distance of the route
    path = convert_arcs_into_path(x_star, n_tot)
    distance = sum(distances[(i,j)]*x_star[i,j] for (i,j) in keys(distances))

    # Output
    path, distance

end



###############
## HEURISTIC ##
###############


function run_heuristic(penalization_factor,
    n, m, r, g, Q, es, ls, I, F,
    neighbours_in, neighbours_out, distances, times,
    display_process)

    # Compute average distance
    avg_distance = 0
    for (i,j) in keys(distances)
        avg_distance += distances[(i,j)]
    end
    avg_distance = avg_distance/m

    # Initialize variables storing results with any initial values
    total_distance = 0
    nb_routes = 0
    routes = []

    # Initialize non-visited set
    non_visited_I = deepcopy(I)
    must_be_visited_I = Set{Int}()

    # Loop while reduced is negative
    if display_process
        while !(isempty(non_visited_I))
            nb_routes += 1
            best_ratio = Inf
            best_distance = Inf
            best_path = []
            println("- ROUTE ", nb_routes)
            println("Generating and finding next best route")
            for i in non_visited_I
                must_be_visited_I = Set{Int}(i)
                path, distance = generate_route(must_be_visited_I,
                    non_visited_I, avg_distance, penalization_factor,
                    n, r, g, Q, es, ls, I, F,
                    neighbours_in, neighbours_out, distances, times)
                nb_newly_visited = length(intersect(non_visited_I, path))
                ratio = distance/nb_newly_visited
                if ratio < best_ratio
                    best_ratio = ratio
                    best_distance = distance
                    best_path = path
                end
            end
            println("Done")
            println("Route: ", best_path)
            println("Route distance: ", best_distance, "\n")
            total_distance += best_distance
            non_visited_I = setdiff(non_visited_I, best_path)
            push!(routes, best_path)
        end
    else
        while !(isempty(non_visited_I))
            nb_routes += 1
            best_ratio = Inf
            best_distance = Inf
            best_path = []
            for i in non_visited_I
                must_be_visited_I = Set{Int}(i)
                path, distance = generate_route(must_be_visited_I,
                    non_visited_I, avg_distance, penalization_factor,
                    n, r, g, Q, es, ls, I, F,
                    neighbours_in, neighbours_out, distances, times)
                nb_newly_visited = length(intersect(non_visited_I, path))
                ratio = distance/nb_newly_visited
                if ratio < best_ratio
                    best_ratio = ratio
                    best_distance = distance
                    best_path = path
                end
            end
            total_distance += best_distance
            non_visited_I = setdiff(non_visited_I, best_path)
            push!(routes, best_path)
        end
    end

    # Output
    total_distance, nb_routes, routes

end
