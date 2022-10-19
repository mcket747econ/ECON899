#So the Helpful functions file has drawshocks, a bellman with interpolation and our parameter structutres
#Need to remember to bring in both it and this model file in our compute code
#Beyond the usual packages, need "Optim" and "Interpolations" otherwise the Bellman won't work


#To Do:
##Full model solver
###Establish structutre
###Figure out how to do the goodness of fit. Seems that we might want an outer while loop.
###Make sure new functions have proper arguments
###Figure out what we want to print. The problem set is pretty useless here, it just says compute an equilibrium
####just print the coefficients?

##Write SimulateCapitalPath
###Policy function wrong right now, need to get the index for k^t_n and Kbar_t and use those in the policy function, not just directly place them
###Make sure all dimensions are correct


## write EstimateRegression
###Returns 4 estimated coefficients
### Need to determine what X and Y are. I think kappa is actually Y but it might be both? MOST IMPORTANT THING TO FIGURE OUT!
###Regression function needs review

##Proper initializer to put everything from the helpful functions together. I just wrote one part of that, I don't know what else needs inititization
using Parameters, Plots, Optim, Interpolations  #last two are needed for the helpful functions
using Random
using Distributions
include("HelpfulFunctions_edits.jl")
# function initialize_mut(R::Results,P::Params) #See PS5 for choices, this could be rewritten as an overall initilize or put into a bigger initilzer
#     @unpack a0, b0, a1, b1, R2 = R
#     @unpack n_k, n_eps, n_K, n_z = P
#     a0 = 0.095
#     b0 = 0.085
#     b1 = 0.999
#     a1 = 0.999
#
#     R2 = [0 0]
#     lambda = 0.1        ##Update with value for lambda
#     pf_k = zeros(n_k, n_eps, n_K,n_z)
#     pf_v = zeros(n_k, n_eps, n_K,n_z)
#     return a0, b0, a1, b1, R2, pf_k, pf_v
# end


function initialize_overall()
    P = Params()
    G = Grids()
    S = Shocks()
    a0 = 0.095
    b0 = 0.085
    b1 = 0.999
    a1 = 0.999

    R2 = [0,0]
    ahat0 = 0
    ahat1 = 0
    bhat0 = 0
    bhat1 = 0
    lambda = 0.1        ##Update with value for lambda
    pf_k = zeros(G.n_k, G.n_eps, G.n_K,G.n_z)
    pf_v = zeros(G.n_k, G.n_eps, G.n_K,G.n_z)
    V_matrix = zeros(P.N,P.T)
    Kappa = zeros(P.T,3)
    R = Results(pf_k,pf_v,a0,a1,b0,b1,R2,ahat0, ahat1, bhat0, bhat1,V_matrix,Kappa)
    return P,G,S,R

end
P,G,S,R = initialize_overall()

# function SimulateCapitalPath(R::Results, P::Params,G::Grids,Epsilon,Z) #need to figure out the proper bounds on this, right now I think its to T-1
#     @unpack pf_k = R
#     @unpack N, T, burn = P
#     @unpack K_grid = G
#
#     kappa::Array{Float64,2} = zeros(T-1,3) #since the first one is Kbar2. I think this is the proper bounds but I am not sure
#     V::Array{Float64,2} = zeros(N,T-1) #Storing the panel. I hate that we have to do this column by column but rows come first
#
#     #do initial column seperately
#     for n in 1:N
#         #Need to get the index for k^t_n and Kbar_t, not directly insert them.
#         K_index = get_index(11.55,K_grid)
#         V[n,1] = pf_k[K_index,Epsilon[n,1],K_index,Z[1]] #11.55 is a hardfix, this probably doesn't work? Need to understand what the policy function is saying better
#     end
#     kappa[1,1] = sum(V[:,1])/N #
#     kappa[1,2] = Z[1]
#     kappa[1,3] = 11.55 #Column three would be capital today and thus X?
#
#     for t in 2:T-1 #I think this is right because in the pseudocode, it only goes to Kbar_T
#         for n in 1:N
#             #See details above about why this doesn't work right
#             V[n,t] = pf_k[V[n,t-1],Epsilon[n,t],kappa[t-1],Z[t]] #need to flip my epsilon I think
#         end
#         kappa[t,1] = sum(V[:,t])/N
#         kappa[t,2] = Z[t] #for the sort later on, tracks the state of the world when the decision was made
#         kappa[t,3] = kappa[t-1,1] #This makes column three yesterday's capital choice, ie today's capital. This is X in that case (I think)
#     end
#
#     #burn the first 1000 columns/rows (kappa is weird)
#     kappa_b::Array{Float64,2} = kappa[burn+1:T,:]
#     V_b::Array{Float64,2} = V[:,burn+1:T]
#
#     return kappa_b, V_b
# end
function SimulateCapitalPath(R::Results, P::Params,G::Grids, Epsilon::Array{Float64, 2},Z::Array{Float64, 1}) #need to figure out the proper bounds on this, right now I think its to T-1
    @unpack pf_k,pf_v = R
    @unpack N, T, burn = P
    @unpack K_grid = G

    kappa::Array{Float64,2} = zeros(T-1,3) #since the first one is Kbar2. I think this is the proper bounds but I am not sure
    V::Array{Float64,2} = zeros(N,T-1) #Storing the panel. I hate that we have to do this column by column but rows come first

    #do initial column seperately
    K_index = trunc(Int64,get_index(11.55,K_grid))
    for n in 1:N
        #Need to get the index for k^t_n and Kbar_t, not directly insert them.
        V[n,1] = pf_k[K_index,trunc(Int64,Epsilon[n,1]),K_index,trunc(Int64,Z[1])] #11.55 is a hardfix, this probably doesn't work? Need to understand what the policy function is saying better
    end
    kappa[1,1] = sum(V[:,1])/N #
    kappa[1,2] = Z[1]
    kappa[1,3] = K_index #Column three would be capital today and thus X?

    for t in 2:T-1 #I think this is right because in the pseudocode, it only goes to Kbar_T
        for n in 1:N
            #See details above about why this doesn't work right
            V[n,t] = pf_k[trunc(Int64,get_index(V[n,t-1],K_grid)),trunc(Int64,Epsilon[n,t]),trunc(Int64,get_index(kappa[t-1,1],K_grid)),trunc(Int64,Z[t])] #need to flip my epsilon I think
        end
        kappa[t,1] = sum(V[:,t])/N
        kappa[t,2] = Z[t] #for the sort later on, tracks the state of the world when the decision was made
        kappa[t,3] = kappa[t-1,1] #This makes column three yesterday's capital choice, ie today's capital. This is X in that case (I think)
    end

    #burn the first 1000 columns/rows (kappa is weird)
    kappa_b::Array{Float64,2} = kappa[burn+1:T-1,:]
    V_b::Array{Float64,2} = V[:,burn+1:T-1]

    return kappa_b, V_b
end

Epsilon, Z = DrawShocks(S,P.N,P.T)
R.pf_k, R.pf_v = Bellman(P,G,S,R)
R.Kappa,R.V_matrix = SimulateCapitalPath(R, P,G, Epsilon,Z)

function state_sort(array) #figure this might be worth doing, will need to check dimensions
    n_b = length(array[:,1])
    good_y::Array{Float64,1} = [0]
    bad_y::Array{Float64,1} = [0]
    good_x::Array{Float64,1} = [0]
    bad_x::Array{Float64,1} = [0]

    for t in 1:n_b #there must be a less janky way to do this but Julia does not like conditional deletes
        if array[t,2] == 1 #assuming that state 1 is good
            good_y = append!(good_y,array[t,1])
            good_x = append!(good_x,array[t,3])
        else
            bad_y = append!(bad_y,array[t,1])
            bad_x = append!(bad_x,array[t,3])
        end
    end
    l_g = length(good_y)
    l_b = length(bad_y)
    good = hcat(good_y[2:l_g],good_x[2:l_g]) #remove the first line and combine the two vectors into a matrix of x and y
    bad = hcat(bad_y[2:l_b],bad_x[2:l_b])

    return good, bad
end

function log_autoreg(x,y,cons::Bool = true) #this one needs checking
    if cons == true
        n_x = length(x)
        lx = log.(x)
        ly = log.(y)
        xi = hcat(ones(n_x),lx) #create a new matrix with a column of ones
        xit = transpose(xi)
        coefs = (xit*xi)^(-1)*(xit*ly) #linear regression
        intercept = coefs[1]
        coef = coefs[2]
        error = ly - xi*coefs
        y_m = ly-ones(n_x)*(ones(n_x)'*ones(n_x))^(-1)*(ones(n_x)'*ly)
        R2 = (transpose(error)*error)/(transpose(y_m)*y_m)
    return intercept, coef, R2
    else
        lx = log.(x)
        ly = log.(y)
        xit = transpose(xi)
        coef = (xit*xi)^(-1)*(xit*ly)
        error = ly - xi*coef
        y_m = ly-ones(n_x)*(ones(n_x)'*ones(n_x))^(-1)*(ones(n_x)'*ly)
        R2 = (transpose(error)*error)/(transpose(y_m)*y_m)
    return coef, R2
    end
end


# function EstimateRegression(kappa,cons::Bool = true )
#     #first do a sort based on z
#     kappa_a, kappa_b = state_sort(kappa) #a is for good state
#
#     #Next run regressions
#     ahat0, ahat1, R2[1] = log_autoreg(kappa_a[:,2],kappa_a[:,1],cons) #not sure what the right Y is, it's not simply epsilon, it might be Kappa again
#     ahat0, ahat1, R2[2] = log_autoreg(kappa_b[:,2],kappa_b[:,1],cons) #not sure what the right Y is
#
#     return ahat0, ahat1, bhat0, bhat1, R2
# end

function EstimateRegression(R::Results,kappa,cons::Bool = true)
    @unpack ahat0, ahat1, bhat0, bhat1, R2 = R
    #first do a sort based on z
    kappa_a, kappa_b = state_sort(kappa) #a is for good state
    #Next run regressions
    ahat0, ahat1, R2[1] = log_autoreg(kappa_a[:,2],kappa_a[:,1],cons) #not sure what the right Y is, it's not simply epsilon, it might be Kappa again
    bhat0, bhat1, R2[2] = log_autoreg(kappa_b[:,2],kappa_b[:,1],cons) #not sure what the right Y is


    return ahat0, ahat1, bhat0, bhat1, R2
end

R.ahat0, R.ahat1, R.bhat0, R.bhat1, R.R2 = EstimateRegression(R,R.Kappa)

function Update_Coef(R::Results,lambda)
    @unpack a0,b0,a1,b1, ahat0,ahat1,bhat0,bhat1 = R
    a0 = lambda .*ahat0 .+ (1 .-lambda) .*a0
    b0 = lambda .*bhat0 .+ (1 .-lambda) .*b0
    a0 = lambda .*ahat1 .+ (1 .-lambda) .*a1
    a0 = lambda .*bhat1 .+ (1 .-lambda) .*b1
    return a0,b0,a1,b1
end

function Solve_KS(P::Params, G::Grids, S::Shocks, R::Results, lambda::Float64=0.5, tol_up = 1e-3, m_gof::Float64 =1 - 1e-2, kill = 10000)
    @unpack a0, b0, a1, b1, R2, pf_k, pf_v = R
    pf_k_temp = zeros(n_k, n_eps, n_K,n_z)
    pf_v_temp = zeros(n_k, n_eps, n_K,n_z)
   #First step is to draw shocks
    Epsilon, Z = DrawShocks(S,P.N,P.T) #use the formula from the helpful functions to get our shocks
    #Second step is to set initial a0,b0,a1,b1

    # R.a0, R.b0, R.a1, R.b1, R.R2, R.pf_k, R.pf_v = initialize_coef(R,P) #I have no idea how initialization works

    stop = 0
    count = 0
    while stop == 0 && minimum(R.R2) < m_gof ##wrong way to do the goodness of fit loop, should be an outer loop but let's get the inner loop working first
        count +=1
        pf_k_temp, pf_pf_v = Bellman(P,G,S,R) #note this returns K' first, different than how we normally do it

        kappa, V_b = SimulateCapitalPath(R,P,Epsilon,Z) #I think this is all correct now
        R.ahat0, R.ahat1, R.bhat0, R.bhat1, R.R2 = EstimateRegression(kappa) #Still need to write this

        if count < kill && (abs(R.ahat0-R.a0)+abs(R.ahat1-R.a1)+abs(R.bhat0-R.b0)+abs(R.bhat1-R.b1)) > tol_up
            R.a0, R.b0, R.a1, R.b1 = Update_Coef(R::Results, lambda) #Still need to write this up
        else
            stop = 1
            println("The coefficients have converged")
            println("Goodness of fit (good, bad): ", R.R2)
        end
    end
    println("Converged in: ",count," iterations")
    println("Predicted a0: ",R.a0)
    println("Predicted a1: ", R.a1)
    println("Predicted b0: ", R.b0)
    println("Predicted b1: ", R.b1)
end

function overall_solve()
    idio_s,agg_s = Drawshocks(S,N,T,sd,)




    stop == 0
    while stop == 0
        pf_k,pf_v = Bellman(P,G,S,R)





    end







end


