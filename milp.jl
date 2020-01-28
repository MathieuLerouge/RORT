using JuMP
using GLPK
m = Model(with_optimizer(GLPK.Optimizer))

## print(m)

## DONNEES ##
n
m
r
g
Q


set_silent(model)
# Z >= 0, PSD
@variable(model, Z[1:m, 1:m], PSD)
@constraint(model, Z .>= 0)
# min Tr(W(I-Z))
@objective(model, Min, tr(W * (Matrix(1.0I, m, m) - Z)))
# Z e = e
@constraint(model, Z * ones(m) .== ones(m))
# Tr(Z) = k
@constraint(model, tr(Z) == k)

JuMP.optimize!(model)
