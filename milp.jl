include("/Users/camillegrange/Documents/RORT/data/data_extraction.jl")

using JuMP
using GLPK
m = Model(with_optimizer(GLPK.Optimizer))


## print(m)

## DONNEES ##

## neighbours_in[i] = liste des j tq existe (j,i)
## neighbours_out[i] = liste des j tq existe (i,j)

n_tot , r, g, Q, es, ls, I, F_prime, neighbours_in, neighbours_out, distances, times = extract_data("/Users/camillegrange/Documents/RORT/data/E_data_1.txt")
n = n_tot - 2

s = [0 for i in 1:n+2]
q = [0 for i in 1:n+2]
C=100000



## VARIABLES ##

@variable(m, x[i in 1:n+2,j in 1:n+2], Bin)
@variable(m, to[i in 1:n+2] >= 0)
@variable(m, u[i in 1:n+2] >= 0)
@variable(m, y[i in 1:n+2] >= 0)

## CONTRAINTES ##


@constraint(m, ct_1[i in I ; i!= 1 && i!=n+2], sum(x[i,j] for j in neighbours_out[i] if j!=1) == 1)
@constraint(m, ct_2[i in F_prime; i!= 1 && i!=n+2], sum(x[i,j] for j in neighbours_out[i] if j!=1) <= 1)
@constraint(m, ct_3[j in 2:n+1], sum(x[j,i] for i in neighbours_out[j] if i != 1) - sum(x[i,j] for i in neighbours_in[j] if i != n+2) == 0)
@constraint(m, ct_4[i in I, j in 2:n+2 ; (i,j) in keys(distances) && i != n+2], to[i]+(times[(i,j)] + s[i])*x[i,j] - ls[1]*(1-x[i,j]) <= to[j])
@constraint(m, ct_5[i in F_prime, j in 2:n+2; (i,j) in keys(distances) && i != 1 && i != n+2], to[i]+times[(i,j)]*x[i,j] + g*(Q-y[i]) - (ls[1] + g*Q)*(1-x[i,j]) <= to[j])
@constraint(m, ct_6[j in 1:n+2], es[j] <= to[j])
@constraint(m, ct_7[j in 1:n+2], to[j] <= ls[j])
@constraint(m, ct_8[i in 1:n+1, j in 2:n+2 ; (i,j) in keys(distances)], u[j] <= u[i] - q[i]*x[i,j] + C*(1-x[i,j]))
@constraint(m, ct_9, u[1] <= C)
@constraint(m, ct_10[i in I, j in 2:n+2; (i,j) in keys(distances) && i != 1 && i != n+2], y[j] <= y[i] - r*distances[(i,j)]*x[i,j] + Q*(1-x[i,j]))
@constraint(m, ct_11[i in F_prime, j in 2:n+2; (i,j) in keys(distances) && i != n+2], y[j] <= Q - r*distances[(i,j)]*x[i,j])


## OBJECTIF ##
@objective(m, Min, sum(distances[(i,j)]*x[i,j] for (i,j) in keys(distances) if i != n+2 && j!= 1))

## RESOLUTION ##

optimize!(m)
println("Objectif : ", objective_value(m))
println("Arcs empruntés :", value.(x))
