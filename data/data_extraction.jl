# Function extract data from data file
# file_path must be a string which contains the extension
function extract_data(file_path)

    # Open the text file and read lines
    f = open(file_path)
    lines = readlines(f)

    # Read first line, i.e general data about the problem
    l = lines[1]
    words = split(l)
    n = parse(Int64, words[1])
    m = parse(Int64, words[2])
    r = parse(Float64, words[3])
    g = parse(Float64, words[4])
    Q = parse(Float64, words[5]) # or Int64?

    # Read first block of lines, i.e data about clients
    es = zeros(Int64, n)
    ls = zeros(Int64, n)
    I = Set{Int}()
    F = Set{Int}()
    for l in lines[2:(n + 1)]
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
    I = setdiff(I, [1, n])
    F = setdiff(F, [1, n])

    # Read second block of lines, i.e data about the arcs
    distances = Dict{Tuple{Int64,Int64},Int64}()
    times = Dict{Tuple{Int64,Int64},Int64}()
    for l in lines[(2+n):end]
        words = split(l)
        i = parse(Int64, words[1])
        j = parse(Int64, words[2])
        dij = parse(Int64, words[3])
        tij = parse(Int64, words[4])
        distances[(i, j)] = dij
        times[(i, j)] = tij
    end

    # Create neighbors sets
    neighbors_in = [Set{Int}() for _ in 1:n]
    neighbors_out = [Set{Int}() for _ in 1:n]
    for (i,j) in keys(distances)
        push!(neighbors_in[j], i)
        push!(neighbors_out[i], j)
    end

    # Close the text file
    close(f)

    # Output all data
    n, m, r, g, Q, es, ls, I, F, neighbors_in, neighbors_out, distances, times

end
