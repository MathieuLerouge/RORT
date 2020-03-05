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

# Let be a vehicule at index i, time t, battery level y.
# The function below compute which customers are feasible from the current state
# where feasible means that the vehicule can reach the customer
# (while satisfying time window, battery and stock)
# and that the vehicule can leave the customer to another node
# (while satisfying the same constraints for this other node)
# Assumption : zero service time, zero demands

function get_feasible_customers(i, t, y, non_visited_customers,
    n_total, r, g, Q, es, ls, stations, neighbors_out, times, distances)

    # Initialize variables
    feasible_customers = Set{Int}()
    second_chance_customers = Set{Int}()

    # Global variable
    l_end_depot = ls[n_total]

    # Loop over non-visited customers j s.t. there is an arc (i,j)
    for j in non_visited_customers
    if (j in neighbors_out[i])

        # Get time and distance to go from i to customer j
        t_ij = times[(i,j)]
        d_ij = distances[(i,j)]

        # Check if reaching the customer j is feasible time-wise
        if (t + t_ij > ls[j])
            # If it is not the case, forget about this customer
            # and study the next one
            continue
        end

        # Check if reaching the customer j is feasible energy-wise
        if (y - r*d_ij < 0)
            # If it is not the case, forget about this customer
            # and study the next one
            # But since if the vehicule would have more battery,
            # the customer would be reachable time-wise,
            # then if the vehicule would go by a recharging station,
            # it may become feasible,
            # therefore we put it in the second-chance customers
            push!(second_chance_customers, j)
            continue
        end

        # Check if leaving the customer j is feasible
        # Initialize a boolean to false
        check_leaving = false

        # Option 1: leaving by going to the end depot
        if (n_total in neighbors_out[j])
            t_j_end = times[(j, n_total)]
            d_j_end = distances[(j, n_total)]
            if (t + t_ij + t_j_end < l_end_depot)
                if (y - r*(d_ij + d_j_end) > 0)
                    check_leaving = true
                else
                    push!(second_chance_customers, j)
                end
            end
        end

        # Option 2: leaving by going by a recharging station
        # Let check this option only if the first one is not possible
        if (!check_leaving)

            # Compute neighbor stations
            neighbor_stations = intersect(stations, neighbors_out[j])

            # If there is no neighboring station
            # then forget about the customer j and study the next one
            if (isempty(neighbor_stations))
                continue
            end

            # Otherwise, let check if there is a feasible_station
            # by looking at each station
            check_feasible_station = false
            while (!isempty(neighbor_stations) && !check_feasible_station)

                # Get one station
                station = pop!(neighbor_stations)

                # Check if the recharging station is reachable
                min_j_station = min(j, station)
                max_j_station = max(j, station)
                t_j_station = times[(min_j_station, max_j_station)]
                d_j_station = distances[(min_j_station, max_j_station)]
                l_station = ls[station]
                if ((t + t_ij + t_j_station > l_station) ||
                    (y - r*(d_ij + d_j_station) < 0))
                        continue
                end

                # Check if leaving the station is feasible
                t_station_end = times[(station, n_total)]
                d_station_end = distances[(station, n_total)]
                arriving_t_station = max(es[station], t + t_ij + t_j_station)
                delta_y_station = Q - (y - r*(d_ij + d_j_station))
                if ((arriving_t_station + g*delta_y_station +
                        t_station_end > l_end_depot) ||
                    (Q - r*d_station_end < 0))
                        continue
                end

                # If both checks
                check_feasible_station = true
            end

            if (!check_feasible_station)
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

# Let be a vehicule at index i, time t, battery level y
# If no non-visited neighboring customers are reachable from i,
# then the vehicule must either go to the end depot or to a recharging station
function choose_end_or_recharge(i, t, y, second_chance_customers,
    n_total, r, g, Q, es, ls, stations, neighbors_out, times, distances)

    # Initialize next node variable
    # next_node will be the index of the end deport or a recharging station
    next_node = -1

    # Case 1: if there is no arc (i,j) with j end depot,
    # then there is no other choice than going to a recharging station
    if (!(n_total in neighbors_out[i]))
        neighbor_stations = intersect(stations, neighbors_out[i])
        if (isempty(neighbor_stations))
            println("Error: no reachable recharging station")
        else
            next_node = pop!(neighbor_stations)
            push!(neighbor_stations, next_node)
        end

    # Case 2: if there exists an arc (i,j) with j end depot,
    # then the vehicule can either go to the end depot
    # or it can go a recharging station in the hope that
    # there will be feasible customers
    else

        # Case 2.1: if there is no feasible customer, even after recharging,
        # which we can know thanks to second_chance_customers,
        # then set next node to end depot
        if (isempty(second_chance_customers))
            next_node = n_total

        # Case 2.2: if there may be feasible customers after recharging
        else

            # Compute neighboring stations
            neighbor_stations = intersect(stations, neighbors_out[i])

            # Case 2.1.1: if there is no neighboring station
            # then set next node to end depot
            if (isempty(neighbor_stations))
                next_node = n_total

            # Case 2.1.2: if there are neighboring stations
            else

                check_feasible_station = false
                while (!isempty(neighbor_stations) && !check_feasible_station)

                    # Get one station
                    station = pop!(neighbor_stations)

                    # Compute leaving time from recharging station
                    arriving_t_station = max(es[station], t)
                    delta_y_station = Q - y
                    leaving_t_station = arriving_t_station + g*delta_y_station

                    # Compute feasible customers from recharging station
                    js = get_feasible_customers(station, leaving_t_station, Q,
                        second_chance_customers,
                        n_total, r, g, Q, es, ls,
                        stations, neighbors_out, times, distances)

                    ## ** WIP!!! **

                    # If there are no feasible customers,
                    # Set next node to end depot
                    if (isempty(js))
                        next_node = n_total

                    # Otherwise,
                    # Set next node to recharging station
                    else
                        next_node = station
                    end
                end
            end
        end
    end

    # Output
    next_node

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
        neighbors_in, neighbors_out, distances, times =
        extract_data(string(data_folder_path, file_name))

    # Initialize sets of non-visited customers and stations
    non_visited_customers = Set(I)
    stations = Set(F)

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
                    non_visited_customers,
                    n_total, r, g, Q, es, ls,
                    stations, neighbors_out, times, distances)

            # Initialize variable for next node
            next_node = -1

            # Case 1: there are feasible customers
            if (!isempty(feasible_customers))
                # Select best customer based on specified criterion
                next_node =
                    select_closest(node, feasible_customers, distances)
                # Remove selected customer from non-visited node
                setdiff!(non_visited_customers, next_node)

            # Case 2 : there is no feasible customers
            else
                # Choose between going to the ending depot
                # or going to a recharging station
                next_node =
                    choose_end_or_recharge(node, time, battery,
                        second_chance_customers,
                        n_total, r, g, Q, es, ls,
                        stations, neighbors_out, times, distances)
            end

            # Update time, battery and total distance
            time = max(time + times[(node, next_node)], es[next_node])
            distance = distances[(node, next_node)]
            total_distance += distance
            battery -= r*distance

            # Update node and data structure
            node = next_node

        end

        # Update variables
        nb_vehicules += 1

    end

    # Output results
    println("Number of vehicules: ", nb_vehicules)
    println("Total distance: ", total_distance)

end




##########
## MAIN ##
##########


# Provide file_name (with extension)
#file_name = "E_data.txt"
file_name = "E_data_1.txt"

# Run greedy heuristic
greedy_heuristic(file_name)
