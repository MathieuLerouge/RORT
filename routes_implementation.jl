# Function which convert representation of the route with binary number coding
# for arcs into a representation of the route as a list of the visited nodes
function convert_arcs_into_path(route, n_tot)
    path = [1]
    i = 1
    while i != n_tot
        j = 2
        while route[i, j] != 1
            j += 1
        end
        i = j
        push!(path, i)
    end
    path
end

function display_route_as_path(route, n_tot)
    println(convert_arcs_into_path(route, n_tot))
end
