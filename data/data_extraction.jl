#####################
## DATA EXTRACTION ##
#####################

## Test

# Function which extracts data from a data file
# file_path must be a string which contains the extension
function extract_data(file_path, nb_duplicates)

    # Open the text file and read lines
    f = open(file_path)
    lines = readlines(f)

    # Read first line, i.e general data about the problem
    l = lines[1]
    words = split(l)
    n_tot = parse(Int64, words[1])
    m = parse(Int64, words[2])
    r = parse(Float64, words[3])
    g = parse(Float64, words[4])
    Q = parse(Float64, words[5]) # or Int64?

    # Read first block of lines, i.e data about clients
    es = zeros(Int64, n_tot)
    ls = zeros(Int64, n_tot)
    I = Set{Int}()
    F = Set{Int}()
    for l in lines[2:(n_tot+1)]
        words = split(l)
        i = parse(Int64, words[1])
        es[i] = parse(Int64, words[2])
        ls[i] = parse(Int64, words[3])
        type = parse(Int64, words[4])
        if (type == 0)
            push!(I, i)
        else
            push!(F, i)
        end
    end
    I = setdiff(I, [1, n_tot])
    F = setdiff(F, [1, n_tot])

    # Read second block of lines, i.e data about the arcs
    distances = Dict{Tuple{Int64,Int64},Int64}()
    times = Dict{Tuple{Int64,Int64},Int64}()
    for l in lines[(2+n_tot):end]
        words = split(l)
        i = parse(Int64, words[1])
        j = parse(Int64, words[2])
        dij = parse(Int64, words[3])
        tij = parse(Int64, words[4])
        distances[(i, j)] = dij
        times[(i, j)] = tij
    end

    # Create neighbors sets
    neighbors_in = [Set{Int}() for _ in 1:n_tot]
    neighbors_out = [Set{Int}() for _ in 1:n_tot]
    for (i,j) in keys(distances)
        push!(neighbors_in[j], i)
        push!(neighbors_out[i], j)
    end

    # Close the text file
    close(f)

    # If we need to dupicate stations
    if nb_duplicates >= 1

        # Put away end depot
        m_prime = m
        end_e = es[n_tot]
        end_l = ls[n_tot]
        end_neighbors_in = neighbors_in[n_tot]
        deleteat!(es, n_tot)
        deleteat!(ls, n_tot)
        deleteat!(neighbors_in, n_tot)
        to_end_distances = Dict{Tuple{Int64,Int64},Int64}()
        to_end_times = Dict{Tuple{Int64,Int64},Int64}()
        for i in end_neighbors_in
            neighbors_out[i] = setdiff(neighbors_out[i], n_tot)
            to_end_distances[(i, n_tot)] = distances[(i, n_tot)]
            to_end_times[(i, n_tot)] = times[(i, n_tot)]
            delete!(distances, (i, n_tot))
            delete!(times, (i, n_tot))
            m_prime -= 1
        end

        # Add duplicates
        n_prime = n_tot-1
        F_prime = copy(F)
        duplicates = Dict{Int64, Set{Int64}}()
        for f in F
            s = Set{Int64}()
            push!(s, f)
            for k in 1:nb_duplicates
                n_prime += 1
                push!(s, n_prime)
                push!(F_prime, n_prime)
                push!(es, es[f])
                push!(ls, ls[f])
                push!(neighbors_in, neighbors_in[f])
                push!(neighbors_out, neighbors_out[f])
                for i in neighbors_in[f]
                    push!(neighbors_out[i], n_prime)
                    distances[(i, n_prime)] = distances[(i, f)]
                    times[(i, n_prime)] = times[(i, f)]
                    m_prime += 1
                end
                for j in neighbors_out[f]
                    push!(neighbors_in[j], n_prime)
                    distances[(n_prime, j)] = distances[(f, j)]
                    times[(n_prime, j)] = times[(f, j)]
                    m_prime += 1
                end
            end
            duplicates[f] = s
        end

        # Put back end depot
        n_prime += 1
        push!(es, end_e)
        push!(ls, end_l)
        push!(neighbors_in, end_neighbors_in)
        for i in end_neighbors_in
            if i in F
                for i_prime in duplicates[i]
                    push!(neighbors_out[i_prime], n_prime)
                    distances[(i_prime, n_prime)] =
                        to_end_distances[(i, n_tot)]
                    times[(i_prime, n_prime)] = to_end_times[(i, n_tot)]
                    m_prime += 1
                end
            else
                push!(neighbors_out[i], n_prime)
                distances[(i, n_prime)] = to_end_distances[(i, n_tot)]
                times[(i, n_prime)] = to_end_times[(i, n_tot)]
                m_prime += 1
            end
        end

        # Update
        n_tot = n_prime
        m = m_prime
        F = F_prime

    end

    # Output all data
    n_tot, m, r, g, Q, es, ls, I, F,
        neighbors_in, neighbors_out, distances, times

end



#####################
## SIMPLIFICATIONS ##
#####################


# Function which simplifies the graph of connections
# by removing unfeasible arcs
function clear_graph(n, m, r, g, Q, es, ls, I, F,
    neighbors_in, neighbors_out, distances, times, s, q, C)

    # Count deleted arcs
    nb_deleted_arcs = 0

    # Delete unfeasible arcs
    for (i,j) in keys(times)

        # Delete arc (i,j) if it violates time windows constraints (14)
        if es[i] + s[i] + times[(i,j)] >= ls[j]
            # Delete arc (i,j) from all data structures
            pop!(times, (i,j))
            pop!(distances, (i,j))
            pop!(neighbors_out[i], j)
            pop!(neighbors_in[j], i)
            nb_deleted_arcs += 1
            continue
        end

        # Delete arc (i,j) if it violates time windows constraints (15)
        if (j,n+2) in keys(times)
            if es[i] + s[i] + times[(i,j)] + s[j] + times[(j,n+2)] >= ls[n+2]
                # Delete arc (i,j) from all data structures
                pop!(times, (i,j))
                pop!(distances, (i,j))
                pop!(neighbors_out[i], j)
                pop!(neighbors_in[j], i)
                nb_deleted_arcs += 1
                continue
            end
        end

        # Delete arc (i,j) if it violates capacity constraints (13)
        if q[i] + q[j] >= C
            # Delete arc (i,j) from all data structures
            pop!(times, (i,j))
            pop!(distances, (i,j))
            pop!(neighbors_out[i], j)
            pop!(neighbors_in[j], i)
            nb_deleted_arcs += 1
            continue
        end

        # Delete arc (i,j) if it violates battery constraints (16)
        should_delete_arc = true
        if (i in I) && (j in I)
            for u in 1:n+1 if (u,i) in keys(times)
                for v in 2:n+2 if (j,v) in keys(times)
                    should_delete_arc = should_delete_arc &&
                        (r*(distances[(u,i)] + distances[(i,j)]
                            + distances[(j,v)]) >= Q)
                end end
            end end
            if should_delete_arc
                # Delete arc (i,j) from all data structures
                pop!(times, (i,j))
                pop!(distances, (i,j))
                pop!(neighbors_out[i], j)
                pop!(neighbors_in[j], i)
                did_delete_arc = true
                nb_deleted_arcs += 1
            end
        end

    end

    # Output
    nb_deleted_arcs

end
