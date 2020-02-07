#############
## INCLUDE ##
#############


# Provide folder_path
#folder_path = "/Users/camillegrange/Documents/RORT/"
folder_path = "/Users/lerougemathieu/Documents/Courses/MPRO/RORT/Project/"

# Do not change
data_folder = "data/"
data_folder_path = string(folder_path, data_folder)
include(string(data_folder_path, "data_extraction.jl"))



#########################
## AUXILIARY FUNCTIONS ##
#########################


    ##################################
    ## FILTERING CUSTOMERS FUNCTION ##
    ##################################

# Assumption : zero service time, zero demands
function get_feasible_customers(i, t, y,
    non_visited_customers, non_visited_stations,
    n_total, r, g, Q, es, ls, neighbors_out, times, distances)

    # Variables
    feasible_customers = Set{Int}()
    second_chance_customers = Set{Int}()
    l_end = ls[n_total]

    # Loop over non-visited customers j s.t. there is an arc (i,j)
    for j in non_visited_customers
    if (j in neighbours_out[i])

        # Get time and distance
        t_ij = times[(i,j)]
        d_ij = distances[(i,j)]

        # Check if reaching the customer j is feasible, time-wise
        if (t + t_ij > ls[j])
            continue
        end

        # Check if reaching the customer j is feasible, energy-wise
        if (y - r*d_ij < 0)
            push!(second_chance_customers, j)
            continue
        end

        # Check if leaving the customer j is feasible
        check_leaving = false

            # Option 1: by going to the end depot
        if (n_total in neighbors_out[j])
            t_jend = times[(j, n_total)]
            d_jend = distances[(j, n_total)]
            if (t + t_ij + t_jend < l_end)
                if (y - r*(d_ij + d_jend) > 0)
                    check_leaving = true
                else
                    push!(second_chance_customers, j)
                end
            end
        end

            # Option 2: by going by a recharging station
        if (!check_leaving)

            # Check if there are non-visited neighbor stations
            neighbor_stations =
                intersect(non_visited_stations, neighbors_out[j])
            if (isempty(neighbor_stations))
                continue
            end

            # If so, get one recharging station
            # By assumption, recharging stations are equivalent
            st = pop!(neighbor_stations)

            # Check if the recharging station is reachable
            min_jst = min(j, st)
            max_jst = max(j, st)
            t_jst = times[(min_jst, max_jst)]
            d_jst = distances[(min_jst, max_jst)]
            l_st = ls[st]
            if ((t + t_ij + t_jst > l_st) || (y - r*(d_ij + d_jst) < 0))
                continue
            end

            # Check if leaving the station is feasible
            t_stend = times[(st, n_total)]
            d_stend = distances[(st, n_total)]
            arriving_t_st = max(es[st], t + t_ij + t_jst)
            delta_y_st = Q - (y - r*(d_ij + d_jst))
            if (arriving_t_st + g*delta_y_st + t_stend > l_end)
                continue
            end
        end

        # Add the customer j to the feasible customer
        push!(feasible_customers, j)

    end
    end

    # Output
    feasible_customers, second_chance_customers

end


    ###################################
    ## ENDING OR RECHARGING FUNCTION ##
    ###################################

function choose_end_or_recharge(i, t, y,
    second_chance_customers, non_visited_stations,
    n_total, r, g, Q, es, ls, neighbors_out, times, distances)

    # Initialize variable
    j_star = -1

    # Case 1: there is not any arc (i,j) with j end depot
    if (!(n_total in neighbors_out[i]))
        neighbor_stations =
            intersect(non_visited_stations, neighbors_out[i])
        if (isempty(neighbor_stations))
            println("Error: no recharging station")
        else
            j_star = pop!(neighbor_stations)
        end

    # Case 2: there exists an arc (i,j) with j end depot
    else

        # Case 2.1 : there are no feasible customers even after recharging
        # Set next node to end depot
        if (isempty(second_chance_customers))
            j_star = n_total

        # Case 2.2 : there may be feasible customers after recharging
        else

            # Check if there are non-visited neighbor stations
            neighbor_stations =
                intersect(non_visited_stations, neighbors_out[i])
            if (isempty(neighbor_stations))
                j_star = n_total

            else
                # If so, get one recharging station
                # By assumption, recharging stations are equivalent
                st = pop!(neighbor_stations)

                # Compute leaving time from recharging station
                arriving_t_st = max(es[st], t)
                delta_y_st = Q - y
                leaving_t_st = arriving_t_st + g*delta_y_st

                # Compute feasible customers from recharging station
                non_visited_stations_st =
                    setdiff(non_visited_stations, st)
                js = get_feasible_customers(st, leaving_t_st, Q,
                    second_chance_customers, non_visited_stations_st,
                    n_total, r, g, Q, es, ls, neighbors_out, times, distances)

                # If there are no feasible customers,
                # Set next node to end depot
                if (isempty(js))
                    j_star = n_total

                # Otherwise,
                # Set next node to recharging station
                else
                    j_star = st
                end
            end
        end
    end

    # Output
    j_star

end


    #########################
    ## SELECTION FUNCTIONS ##
    #########################

function select_closest(i, feasible_customers, distances)

    # Initialize variables
    j_star = -1
    d_star_ij = Inf

    # Find the closest customer
    for j in feasible_customers
        d_ij = distances[(i,j)]
        if (d_ij < d_star_ij)
            j_star = j
            d_star_ij = d_ij
        end
    end

    # Output
    j_star

end


# ** To do **
# If positive demands
function select_largest_ratio_distance_over_demand(
    i, feasible_customers, distances)

    # Initialiaze variables
    j_star = -1
    d_star_ij = Inf

    # Find the largest ratio distance over demand
    # ** To do **

    # Output
    j_star

end




######################
## GREEDY HEURISTIC ##
######################


# Assumption : zero service time, zero demands
function greedy_heuristic(file_name)

    # Extract data
    n_total, m, r, g, Q, es, ls, I, F,
        neighbours_in, neighbors_out, distances, times =
        extract_data(string(data_folder_path, file_name))

    # Initialize data structures
    non_visited_customers = Set(I)
    non_visited_stations = Set(F)

    # Initialize variables
    nb_vehicules = 0
    total_distance = 0

    # While all customers have not been visited
    while (!isempty(non_visited_customers))

        # Initialiaze route
        node = 1
        time = es[node]
        battery = Q

        # While vehicule is not at the end depot
        while (node != n_total)

            # Compute feasible customers
            feasible_customers, second_chance_customers =
                get_feasible_customers(node, time, battery,
                non_visited_customers, non_visited_stations,
                n_total, r, g, Q, es, ls, neighbors_out, times, distances)

            # Initialize case variable
            case = 3

            # Case 1: there are feasible customers
            if (!isempty(feasible_customers))
                # Select best customer based on specified criterion
                next_node =
                    select_closest(node, feasible_customers, distances)
                # Set case
                case = 1

            # Cases 2 or 3 : choose between ending route or recharging
            else
                next_node =
                    choose_end_or_recharge(node, time, battery,
                    second_chance_customers, non_visited_stations,
                    n_total, r, g, Q, es, ls, neighbors_out, times, distances)

                # Case 2: going to recharging station
                if (next_node != n_total)
                    # Set case
                    case = 2

                # Case 3: going to end depot
                else
                    # Set case
                    case = 3
                end
            end

            # Update time, battery and total distance
            time = max(time + times[(node, next_node)], es[next_node])
            distance = distances[(node, next_node)]
            total_distance += distance
            battery -= r*distance

            # Update node and data structure
            node = next_node
            if (case == 1)
                setdiff!(non_visited_customers, node)
            elseif (case == 2)
                setdiff!(non_visited_stations, node)
            end

            # Temporary
            #node = n_total
        end

        # Update variables
        nb_vehicules += 1

        # Temporary
        empty!(non_visited_customers)

    end

    # Output results
    println("Number of vehicules: ", nb_vehicules)
    println("Total distance: ", total_distance)

end




##########
## MAIN ##
##########


# Provide file_name (with extension)
file_name = "E_data.txt"
#file_name = "E_data_1.txt"

# Run greedy heuristic
greedy_heuristic(file_name)
