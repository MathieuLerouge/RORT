################
## LIBRAIRIES ##
################


using DataStructures



###############
## STRUCTURE ##
###############


# Set global variables
start_depot_index = 1
tolerance = 0.000001

# Define a structure
struct Tag
    time::Float64
    distance::Float64
    battery::Float64
    predecessor::Int
end



###############
## FUNCTIONS ##
###############


# Function to initialize the tags stored at the nodes
function initialize_nodes_tags(n_tot, r, g, Q, es, I, F)

    # Initialize list of tags
    tags = [[Tag(es[start_depot_index], 0, Q, 0)]]
    deltas = [[Tag(es[start_depot_index], 0, Q, 0)]]
    for i in 2:n_tot
        push!(tags, [])
        push!(deltas, [])
    end
    for i in 2:n_tot
        push!(tags[i], Tag(ls[i], Inf, 0, 0))
    end

    # Output
    tags, deltas

end


# Function to initializa the queue of nodes that must be explored
function initialize_nodes_exploration_queue()

    # Initialize exploration queue
    exploration_queue = Queue{Int}()
    enqueue!(exploration_queue, start_depot_index)

    # Output
    exploration_queue

end


# Function to compute the sequential update of the tags tags_j
# of a successor j of a node i based on a translation of the tags delta_i
function compute_sequential_update(i, i_is_a_station, delta_i, tags_j, e_j, l_j,
    d_ij, t_ij, c_ij, r, g, Q)

    # Initialize delta_j
    delta_j = []

    # Study each potential new tag obtained by translating a tag of delta_i
    for tag_i in delta_i

        # Compute the new tag obtained by translation
        d_i = tag_i.distance
        t_i = tag_i.time
        y_i = tag_i.battery
        if (i_is_a_station)
            t_i += g*(Q - y_i)
            y_i = Q
        end
        y_j = y_i - r*c_ij
        t_j = max(e_j, t_i + t_ij)
        d_j = d_i + d_ij
        new_tag_j = Tag(t_j, d_j, y_j, i)

        # Check if the new time t_j fits in the time window
        # and if the charge y_j is positive
        if (t_j <= l_j) #&& (y_j >= 0) # ** WIP **

            # Set variables
            must_insert_tag = false
            nb_tags_j = length(tags_j)

            # Find the tag k, named key tag, of the successor j which has the
            # smallest time t_j[k] that is strictly superior to the new time t_j
            # if such k exists
            # (in general case, 1 < k <= nb_tags and
            # k is such that t_j[k-1] <= t_j < t_j[k])
            k = 1
            while (k <= nb_tags_j) && (t_j >= tags_j[k].time)
                k += 1
            end

            # Case 0:
            # When 1 < k, new time t_j is equal to the time of tag k-1
            # (i.e. when 1 < k <= nb_tags, t_j[k-1] = t_j < t_j[k])
            if (k > 1) && (tags_j[k-1].time == t_j)

                # Case 0.a:
                # New distance d_j is strictly smaller than
                # the distance d_j[k-1] of tag k-1
                # (i.e. t_j[k-1] = t_j and d_j[k-1] > d_j)
                if (tags_j[k-1].distance > d_j)

                    # Delete tag k-1
                    deleteat!(tags_j, k-1)

                    # Update
                    k -= 1
                    nb_tags_j -= 1

                    # Flag the need to insert the new tag
                    must_insert_tag = true

                # Case b:
                # New distance d_j is larger or equal to
                # the distance d_j[k-1] of tag k-1
                # (i.e. t_j[k-1] = t_j and d_j[k-1] <= d_j)
                else

                    # Skip the new tag
                    continue

                end
            end

            # From here, 1 <= k <= nb_tags + 1
            # and if 1 < k <= nb_tags then t_j[k-1] < t_j < t_j[k]

            # Case 1:
            # There is no tag at node j
            if (nb_tags_j == 0)

                # Add the new tag
                push!(tags_j, new_tag_j)

                # Save the change
                push!(delta_j, new_tag_j)

                # Flag no need to insert the new tag
                must_insert_tag = false

            # Case 2:
            # New time t_j is strictly larger than any other time but
            # (i.e. t_j[-1] < t_j)
            elseif (nb_tags_j < k)

                # Case 2.a:
                # New distance d_j is strictly smaller than
                # the distance of the last tag
                # (i.e. t_j[-1] < t_j and d_j[-1] > d_j)
                if (tags_j[nb_tags_j].distance > d_j)

                    # Add the new tag
                    push!(tags_j, new_tag_j)

                    # Save the change
                    push!(delta_j, new_tag_j)

                    # Flag no need to insert the new tag
                    must_insert_tag = false

                # Case 2.b:
                # New distance d_j is larger or equal to
                # the distance of the last tag
                # (i.e. t_j[-1] < t_j and d_j[-1] <= d_j)
                else

                    # Skip the new tag
                    continue

                end

            # From here, 1 <= k <= nb_tags
            # - if k = 1 then t_j < t_j[1]
            # - if 1 < k <= nb_tags then t_j[k-1] < t_j < t_j[k]

            # Case 3:
            # New distance d_j is smaller than the one of the key tag d_j[k]
            # (i.e. t_j[k-1] < t_j < t_j[k] and d_j <= d_j[k])
            elseif (d_j <= tags_j[k].distance)

                # Remove any tag which has a larger time and a larger distance
                while (k <= nb_tags_j) && (d_j <= tags_j[k].distance)
                    deleteat!(tags_j, k)
                    nb_tags_j -= 1
                end

                # Border effect
                if (nb_tags_j == 0)

                    # Add the new tag
                    push!(tags_j, new_tag_j)

                    # Save the change
                    push!(delta_j, new_tag_j)

                    # Flag no need to insert the new tag
                    must_insert_tag = false

                else

                    # Flag the need to insert the new tag
                    must_insert_tag = true

                end

            # From here, 1 <= k <= nb_tags
            # - if k = 1 then t_j < t_j[1] and d_j > d_j[1]
            # - if 1 < k <= nb_tags then t_j[k-1] < t_j < t_j[k] and d_j > d_j[k]

            # Case 4:
            # New distance d_j is not smaller than the one of the key tag d_j[k]
            # but is smaller than the one of the tag before
            # (i.e. t_j[k-1] < t_j < t_j[k] and d_j[k-1] > d_j > d_j[k])
            elseif (k == 1) || ((k > 1) && (tags_j[k-1].distance) > d_j)

                # Flag the need to insert the new tag
                must_insert_tag = true

            end

            # If the need to insert has been raised
            if must_insert_tag

                # Insert tag at the right location
                insert!(tags_j, k, new_tag_j)

                # Save the change
                push!(delta_j, new_tag_j)

            end
        end
    end

    # ** WIP **

    # Output updated tags of the successor j
    tags_j, delta_j

end


# Function to compute the tags of the successors js of a node i
function compute_successors_new_tags(index_i, i_is_a_station, delta_i,
    successors_js, successors_tags, successors_es, successors_ls,
    successors_distances, successors_times, successors_costs, r, g, Q)

    # Initialize data structure about modified successors
    successors_js = collect(successors_js)
    modified_successors_js = []
    modified_successors_tags = []
    modified_successors_incr_deltas = []

    # Update each successor's tags
    for j in 1:length(successors_js)

        # Compute updated tags of successor j
        tags_j, delta_j = compute_sequential_update(index_i, i_is_a_station,
            delta_i, successors_tags[j], successors_es[j], successors_ls[j],
            successors_distances[j], successors_times[j], successors_costs[j],
            r, g, Q)

        # If there are changes due to the update
        # then save the changes
        if (length(delta_j) > 0)
            push!(modified_successors_js, successors_js[j])
            push!(modified_successors_tags, tags_j)
            push!(modified_successors_incr_deltas, delta_j)
        end

    end

    # Output
    modified_successors_js, modified_successors_tags,
        modified_successors_incr_deltas

end


function add_to_nodes_exploration_queue(nodes_exploration_queue,
    modified_successors_js, sorting_keys)

    # Sort modified successors by keys
    sorted_pairs = sort(collect(zip(sorting_keys, modified_successors_js)))

    # Add the successors to the exploration queue
    for pair in sorted_pairs
        enqueue!(nodes_exploration_queue, pair[2])
    end

    # Output
    nodes_exploration_queue

end


# Function which finds the right predecessor tag
function get_predecessor_tag(predecessor_distance, predecessor_time,
    predecessor_tags)

    # Find right predecessor tag
    right_predecessor_tag_index = 1
    did_find_right_tag = false
    while (right_predecessor_tag_index <= length(predecessor_tags)) &&
        (!did_find_right_tag)
            predecessor_tag = predecessor_tags[right_predecessor_tag_index]
            if (predecessor_tag.time <= predecessor_time + tolerance) &&
                (abs(predecessor_tag.distance - predecessor_distance) < tolerance)
                did_find_right_tag = true
            else
                right_predecessor_tag_index += 1
            end
    end
    predecessor_tag = predecessor_tags[right_predecessor_tag_index]

    # Output
    predecessor_tag

end


# Function to run the algorithm of shortest path with time windows provided
# in the article
function run_constrained_shortest_path(n, r, g, Q, es, ls, I, F,
    neighbours_out, distances, times, costs)

    # Define variable
    n_tot = n+2

    # Initialize variables
    nodes_tags, nodes_deltas = initialize_nodes_tags(n_tot, r, g, Q, es, I, F)
    nodes_exploration_queue = initialize_nodes_exploration_queue()

    # Explore and update tags
    while !(isempty(nodes_exploration_queue))

        # Select a node and get data about its successors
        current_node = dequeue!(nodes_exploration_queue)
        current_node_is_a_station = (current_node in F)
        successors_js = neighbours_out[current_node]
        successors_tags = []
        successors_es = []
        successors_ls = []
        successors_distances = []
        successors_times = []
        successors_costs = []
        for j in successors_js
            push!(successors_tags, nodes_tags[j])
            push!(successors_es, es[j])
            push!(successors_ls, ls[j])
            push!(successors_distances, distances[(current_node,j)])
            push!(successors_times, times[(current_node,j)])
            push!(successors_costs, costs[(current_node,j)])
        end

        # Update tags and deltas of modified successors
        modified_successors_js, modified_successors_tags,
            modified_successors_incr_deltas =
            compute_successors_new_tags(current_node, current_node_is_a_station,
                nodes_deltas[current_node], successors_js, successors_tags,
                successors_es, successors_ls, successors_distances,
                successors_times, successors_costs, r, g, Q)
        sorting_keys = []
        counter = 1
        for j in modified_successors_js
            nodes_tags[j] = modified_successors_tags[counter]
            append!(nodes_deltas[j], modified_successors_incr_deltas[counter])
            push!(sorting_keys, es[j])
            counter += 1
        end

        # Add modified successors to the exploration queue
        nodes_exploration_queue =
            add_to_nodes_exploration_queue(nodes_exploration_queue,
                modified_successors_js, sorting_keys)

        # Empty delta of current node
        nodes_deltas[current_node] = []

    end

    # Get shortest distance
    shortest_distance = Inf
    best_tag_index = 0
    for k in 1:length(nodes_tags[n_tot])
        tag = nodes_tags[n_tot][k]
        if (tag.distance < shortest_distance)
            shortest_distance = tag.distance
            best_tag_index = k
        end
    end

    # Get shortest path
    shortest_path = [n_tot]
    tag = nodes_tags[n_tot][best_tag_index]
    node = n_tot
    lattest_time = ls[n_tot]
    predecessor = tag.predecessor
    while (predecessor != start_depot_index)
        predecessor_distance = tag.distance - distances[(predecessor, node)]
        predecessor_time =
            min(lattest_time - times[(predecessor, node)], ls[predecessor])
        tag = get_predecessor_tag(predecessor_distance, predecessor_time,
            nodes_tags[predecessor])
        node = predecessor
        predecessor = tag.predecessor
        pushfirst!(shortest_path, node)
    end
    pushfirst!(shortest_path, predecessor)

    # Output
    shortest_path, shortest_distance

end
