################
## LIBRAIRIES ##
################

using JuMP
using GLPK



##########
## MILP ##
##########


function run_MILP(n, r, g, Q, es, ls, I, F_prime,
    neighbours_in, neighbours_out, distances, times, s, q, C)

    ## MODEL ##

    # Define the model
    model = Model(GLPK.Optimizer)


    ## DECISION VARIABLES ##

    # Define the decision variables
    @variable(model, x[i in 1:n+1,j in 2:n+2], Bin)
    @variable(model, tau[i in 1:n+2] >= 0)
    # Assumption: no demand to deliver
    #@variable(model, u[i in 1:n+2] >= 0)
    @variable(model, y[i in 1:n+2] >= 0)


    ## CONSTRAINTS ##

    # Define the constraints

    # - Delivery constraint (3)
    @constraint(model, ct_2[i in I],
        sum(x[i,j] for j in neighbours_out[i] if j!=1) == 1)

    # - Recharging station constraint (3)
    @constraint(model, ct_3[i in F_prime],
        sum(x[i,j] for j in neighbours_out[i] if j!=1) <= 1)

    # - Flow conservation (4)
    @constraint(model, ct_4[j in 2:n+1],
        sum(x[j,i] for i in neighbours_out[j] if i != 1) -
        sum(x[i,j] for i in neighbours_in[j] if i != n+2) == 0)

    # - Time feasability when leaving customers and depot (5)
    @constraint(model, ct_5[i in union(I, 1), j in 2:n+2 ; (i,j) in keys(distances)],
        tau[i]+(times[(i,j)] + s[i])*x[i,j] - ls[1]*(1-x[i,j]) <= tau[j])

    # - Time feasability when leaving recharging stations (6)
    if (!isempty(F_prime))
        @constraint(model, ct_6[i in F_prime, j in 2:n+2;
            (i,j) in keys(distances)],
            tau[i]+times[(i,j)]*x[i,j] + g*(Q-y[i]) - (ls[1] + g*Q)*(1-x[i,j])
            <= tau[j])
    end

    # - Time windows (7)
    @constraint(model, ct_7[j in 1:n+2], es[j] <= tau[j])
    @constraint(model, ct_7bis[j in 1:n+2], tau[j] <= ls[j])

    # - Demand feasability (8) and (9)
    # Assumption: no demand to deliver
    #@constraint(model, ct_8[i in 1:n+1, j in 2:n+2 ; (i,j) in keys(distances)],
    #    u[j] <= u[i] - q[i]*x[i,j] + C*(1-x[i,j]))
    #@constraint(model, ct_9, u[1] <= C)

    # - Battery feasability (10) and (11)
    @constraint(model, ct_10[i in I, j in 2:n+2;
        (i,j) in keys(distances)],
        y[j] <= y[i] - r*distances[(i,j)]*x[i,j] + Q*(1-x[i,j]))
    @constraint(model, ct_11[i in union(F_prime,1), j in 2:n+2;
        (i,j) in keys(distances) && i != n+2],
        y[j] <= Q - r*distances[(i,j)]*x[i,j])

    # Define the objective function
    @objective(model, Min,
        sum( distances[(i,j)]*x[i,j]
            for (i,j) in keys(distances) if i != n+2 && j!= 1))


    ## SOLVING PROCESS ##

    # Run the solving process
    solving_time = @elapsed optimize!(model)


    ## RESULTS ##

    # Get solutions
    xs_star = value.(x)
    total_distance = objective_value(model)

    # Compute number of routes/vehicules
    nb_vehicules = Int(sum(xs_star[1,:]))

    # Outputs
    total_distance, nb_vehicules, solving_time

end
