## LIBRAIRIES ##

using JuMP
using GLPK


## INCLUDE ##

#folder_path = "/Users/camillegrange/Documents/RORT/"
folder_path = "/Users/lerougemathieu/Documents/Courses/MPRO/RORT/Project/"
data_folder = "data/"
data_folder_path = string(folder_path, data_folder)
include(string(data_folder_path, "data_extraction.jl"))


## MODEL ##

model = Model(with_optimizer(GLPK.Optimizer))


## DATA ##

## neighbours_in[i] = liste des j tq existe (j,i)
## neighbours_out[i] = liste des j tq existe (i,j)

file_name = "E_data.txt"
#file_name = "E_data_1.txt"
n_tot, m, r, g, Q, es, ls, I, F_prime,
    neighbours_in, neighbours_out, distances, times =
    extract_data(string(data_folder_path, file_name))

n = n_tot - 2

s = [0 for i in 1:n+2]
q = [0 for i in 1:n+2]
C = 100000


## VARIABLES ##

@variable(model, x[i in 1:n+1,j in 2:n+2], Bin)
@variable(model, tau[i in 1:n+2] >= 0)
@variable(model, u[i in 1:n+2] >= 0)
@variable(model, y[i in 1:n+2] >= 0)


## CONSTRAINTS ##

@constraint(model, ct_1[i in I ; i!= 1 && i!=n+2],
    sum(x[i,j] for j in neighbours_out[i] if j!=1) == 1)
@constraint(model, ct_2[i in F_prime; i!= 1 && i!=n+2],
    sum(x[i,j] for j in neighbours_out[i] if j!=1) <= 1)
@constraint(model, ct_3[j in 2:n+1],
    sum(x[j,i] for i in neighbours_out[j] if i != 1) -
    sum(x[i,j] for i in neighbours_in[j] if i != n+2) == 0)
@constraint(model, ct_4[i in I, j in 2:n+2 ; (i,j) in keys(distances) && i != n+2],
    tau[i]+(times[(i,j)] + s[i])*x[i,j] - ls[1]*(1-x[i,j]) <= tau[j])
@constraint(model, ct_5[i in F_prime, j in 2:n+2;
    (i,j) in keys(distances) && i != 1 && i != n+2],
    tau[i]+times[(i,j)]*x[i,j] + g*(Q-y[i]) - (ls[1] + g*Q)*(1-x[i,j]) <= tau[j])
@constraint(model, ct_6[j in 1:n+2], es[j] <= tau[j])
@constraint(model, ct_7[j in 1:n+2], tau[j] <= ls[j])
@constraint(model, ct_8[i in 1:n+1, j in 2:n+2 ; (i,j) in keys(distances)],
    u[j] <= u[i] - q[i]*x[i,j] + C*(1-x[i,j]))
@constraint(model, ct_9, u[1] <= C)
@constraint(model, ct_10[i in I, j in 2:n+2;
    (i,j) in keys(distances) && i != 1 && i != n+2],
    y[j] <= y[i] - r*distances[(i,j)]*x[i,j] + Q*(1-x[i,j]))
@constraint(model, ct_11[i in F_prime, j in 2:n+2;
    (i,j) in keys(distances) && i != n+2],
    y[j] <= Q - r*distances[(i,j)]*x[i,j])


## OBJECTIF ##

@objective(model, Min,
    sum(distances[(i,j)]*x[i,j] for (i,j) in keys(distances) if i != n+2 && j!= 1))


## SOLVING ##

optimize!(model)


## RESULTS ##

xs_star = value.(x)
nb_vehicules = Int(sum(xs_star[1,:]))
println("Number of vehicules: ", nb_vehicules)
println("Total distance: ", objective_value(model))
#println("Arcs:", xs_star)
