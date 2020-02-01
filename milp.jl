include("/Users/camillegrange/Documents/RORT/data/data_extraction.jl")

using JuMP
using GLPK
m = Model(with_optimizer(GLPK.Optimizer))


## print(m)

## DONNEES ##

#n , r , g ,Q  = extract_data("/Users/camillegrange/Documents/RORT/data/E_data.txt")
## I, F_prime, distances, times, es, ls, s, q, neighbours_in, neighbours_out
## neighbours_in[i] = liste des j tq existe (j,i)
## neighbours_out[i] = liste des j tq existe (i,j)

n = 3
r = 1
g = 1
Q = 100000000
I = [2, 3, 4]
F_prime = []
distances = Dict()
distances[(1,2)] = 1
distances[(2,3)] = 1
distances[(1,3)] = 1
distances[(2,5)] = 1
distances[(3,5)] = 1
distances[(4,5)] = 1
distances[(1,4)] = 1
times = Dict()
times[(1,2)] = 1
times[(2,3)] = 1
times[(1,3)] = 1
times[(2,5)] = 1
times[(3,5)] = 1
times[(4,5)] = 1
times[(1,4)] = 1
es = [0 for i in 1:5]
ls = [10 for i in 1:5]
s = [0 for i in 1:5]
q = [0 for i in 1:5]
neighbours_out = [[2,3,4],[3,5],[5],[5],[]]
neighbours_in = [[],[1],[2,1],[1],[2,3,4]]
C = 1000


## VARIABLES ##

@variable(m, x[i in 1:n+2,j in 1:n+2], Bin)
@variable(m, to[i in 1:n+2] >= 0)
@variable(m, u[i in 1:n+2] >= 0)
@variable(m, y[i in 1:n+2] >= 0)

## CONTRAINTES ##

@constraint(m, ct_1[i in I], sum(x[i,j] for j in neighbours_out[i]) == 1)
@constraint(m, ct_2[i in F_prime], sum(x[i,j] for j in neighbours_out[i]) <= 1)
@constraint(m, ct_3[j in 2:n+1], sum(x[j,i] for i in neighbours_out[j]) - sum(x[i,j] for i in neighbours_in[j]) == 0)
@constraint(m, ct_4[i in I, j in 2:n+2 ; (i,j) in keys(distances)], to[i]+(times[(i,j)] + s[i])*x[i,j] - ls[1]*(1-x[i,j]) <= to[j])
@constraint(m, ct_4_bis[i in 1:1, j in 2:n+2 ; (i,j) in keys(distances)], to[i]+(times[(i,j)] + s[i])*x[i,j] - ls[1]*(1-x[i,j]) <= to[j])

@constraint(m, ct_5[i in F_prime, j in 2:n+2, (i,j) in keys(distances)], to[i]+times[i,j]*x[i,j] + g*(Q-y[i]) - (ls[1] + g*Q)*(1-x[i,j]) <= to[j])
@constraint(m, ct_6[j in 1:n+2], es[j] <= to[j])
@constraint(m, ct_7[j in 1:n+2], to[j] <= ls[j])
#@constraint(m, ct_8[i in 1:n+1, j in 2:n+2 ; (i,j) in keys(distances)], u[j] <= u[i] - q[i]*x[i,j] + C(1-x[i,j]))
@constraint(m, ct_9, u[1] <= C)
@constraint(m, ct_10[i in I, j in 2:n+2; (i,j) in keys(distances)], y[j] <= y[i] - r*distances[(i,j)]*x[i,j] + Q*(1-x[i,j]))
@constraint(m, ct_11[i in F_prime, j in 2:n+2; (i,j) in keys(distances)], y[j] <= Q - r*distances[(i,j)]*x[i,j])
@constraint(m, ct_11_bis[i in 1:1, j in 2:n+2; (i,j) in keys(distances)], y[j] <= Q - r*distances[(i,j)]*x[i,j])


## OBJECTIF ##
@objective(m, Min, sum(distances[(i,j)]*x[i,j] for (i,j) in keys(distances)))

## RESOLUTION ##

optimize!(m)
println("Objectif : ", objective_value(m))
println("Arcs empruntés :", value.(x))
