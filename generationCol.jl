## LIBRAIRIES ##

using JuMP
using GLPK


## INCLUDE ##

folder_path = "/Users/camillegrange/Documents/RORT/"
#folder_path = "/Users/lerougemathieu/Documents/Courses/MPRO/RORT/Project/"
data_folder = "data/"
data_folder_path = string(folder_path, data_folder)
include(string(data_folder_path, "data_extraction.jl"))


## DATA ##

## neighbours_in[i] = liste des j tq existe (j,i)
## neighbours_out[i] = liste des j tq existe (i,j)

file_name = "E_data.txt"
#file_name = "E_data_1.txt"
n_tot, m, r, g, Q, es, ls, I, F_prime,
    neighbours_in, neighbours_out, distances, times =
    extract_data(string(data_folder_path, file_name))

n = n_tot - 2

#Pas de temps de services :
s = [0 for i in 1:n+2]
#Pas de colis à livrer :
q = [0 for i in 1:n+2]
#Capacité des voitures :
C=100000

function solve_slave(var_duales)

    model = Model(with_optimizer(GLPK.Optimizer))


    ## VARIABLES ##

    @variable(model, x[i in 1:n+1,j in 2:n+2], Bin)
    @variable(model, tau[i in 1:n+2] >= 0)
    @variable(model, u[i in 1:n+2] >= 0)
    @variable(model, y[i in 1:n+2] >= 0)


    ## CONSTRAINTS ##

    @constraint(model, ct_2[i in F_prime],
        sum(x[i,j] for j in neighbours_out[i] if j!=1) <= 1)
    @constraint(model, ct_3[j in 2:n+1],
        sum(x[j,i] for i in neighbours_out[j] if i != 1) -
        sum(x[i,j] for i in neighbours_in[j] if i != n+2) == 0)
    @constraint(model, ct_4[i in union(I, 1), j in 2:n+2 ; (i,j) in keys(distances)],
        tau[i]+(times[(i,j)] + s[i])*x[i,j] - ls[1]*(1-x[i,j]) <= tau[j])
    @constraint(model, ct_5[i in F_prime, j in 2:n+2;
        (i,j) in keys(distances)],
        tau[i]+times[(i,j)]*x[i,j] + g*(Q-y[i]) - (ls[1] + g*Q)*(1-x[i,j]) <= tau[j])
    @constraint(model, ct_6[j in 1:n+2], es[j] <= tau[j])
    @constraint(model, ct_7[j in 1:n+2], tau[j] <= ls[j])
    @constraint(model, ct_8[i in 1:n+1, j in 2:n+2 ; (i,j) in keys(distances)],
        u[j] <= u[i] - q[i]*x[i,j] + C*(1-x[i,j]))
    @constraint(model, ct_9, u[1] <= C)
    @constraint(model, ct_10[i in I, j in 2:n+2;
        (i,j) in keys(distances)],
        y[j] <= y[i] - r*distances[(i,j)]*x[i,j] + Q*(1-x[i,j]))
    @constraint(model, ct_11[i in union(F_prime,1), j in 2:n+2;
        (i,j) in keys(distances) && i != n+2],
        y[j] <= Q - r*distances[(i,j)]*x[i,j])


    ## OBJECTIF ##

    @objective(model, Min,
        sum(distances[(i,j)]*x[i,j] for (i,j) in keys(distances)) - sum(var_duales[i]*x[i,j] for i in I for j in 2:n+2 if (i,j) in keys(distances)))

    ## SOLVING ##

    optimize!(model)


    ## RESULTS ##

    xs_star = value.(x)

    nb_vehicules = Int(sum(xs_star[1,:]))
    println("Number of vehicules: ", nb_vehicules)
    println("Total distance: ", objective_value(model))
    #println("Arcs:", xs_star)
    xs_star, x
end

function solve_master(routes)

    model = Model(with_optimizer(GLPK.Optimizer))
    K = length(routes)

    ## VARIABLES ##
    @variable(model, delta[i in 1:K], Bin)

    ## CONSTRAINTS ##

    @constraint(model, ct_1[i in I],
        sum(sum(routes[k][i,j] for j in 2:n+2 if (i,j) in keys(distances))*delta[k] for k in 1:K) == 1)

    ## OBJECTIF ##

    @objective(model, Min,
        sum(sum(routes[k][i,j]*distances[(i,j)] for i in 1:n+1 for j in 2:n+2 if (i,j) in keys(distances))*delta[k] for k in 1:K))

    ## SOLVING ##

    optimize!(model)


    ## RESULTS ##

    obj = value.(delta)
    dual_values = getRowDualGLPK(model)
    nb_routes = Int(sum(delta_star[1,:]))
    println("Number of routes: ", nb_routes)
    println("Total distance: ", objective_value(model))
    #println("Routes:", delta_star)
    obj, delta, dual_values
end

function fill_graph(n_tot, m, r, g, Q, es, ls, I, F_prime,
    neighbours_in, neighbours_out, distances, times)
    infinite = 10000
    for i in I
        neighbours_in[i] = union(neighbours_in[i], 1)
        neighbours_out[i] = union(neighbours_out[i], n+2)
        if((1,i) in keys(distances))
            continue
        else
            distances[(1,i)] = infinite
            times[(1,i)] = infinite
        end
        if((i, n+2) in keys(distances))
            continue
        else
            distances[(i, n+2)] = infinite
            times[(i,n+2)] = infinite
        end
    end
end

function routes_init()
    khi = []
    for i in I
        route_i = zeros(Int64, n+2, n+2)
        route_i[1,i] = 1 # l'arc (i,j) est à l'indice (i-1)*(n+2) + j (lignes commencent à 0 mais colonnes à 1)
        route_i[i,n+2] = 1
        append!(khi, route_i)
    end
    khi
end

function genCol(routes_ini)

    redCost = -1
    routes = deepcopy(routes_ini)

    while(redCost < 0)

        obj, delta, var_duales = solve_master(routes)
        redCost, newRoute = solve_slave(var_duales)
        append!(routes, newRoute)
    end
    routes, obj, delta
end


#####################"
#       MAIN        #
###################"#

fill_graph(n_tot, m, r, g, Q, es, ls, I, F_prime, neighbours_in, neighbours_out, distances, times)

routes_ini = routes_init()

genCol(routes_ini)
