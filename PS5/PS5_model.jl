#So the Helpful functions file has drawshocks, a bellman with interpolation and our parameter structutres
#Need to remember to bring in both it and this model file in our compute code
#Beyond the usual packages, need "Optim" and "Interpolations" otherwise the Bellman won't work 


#To Do:  
##Full model solver
###Establish structutre
###Put in precomposed functions (mostly done)
###Figure out how to do the goodness of fit. Seems that we might want an outer while loop. 
###Make sure new functions have proper arguments

##Write SimulateCapitalPath
###This needs to be run now, but should work fine? 
###Make sure all dimensions are correct 
###Add in another column to Kappa so I track K today vs K tomorrow right? For regression 


## write EstimateRegression
###Returns 4 estimated coefficients 
### Need to determine what X and Y are. I think kappa is actually Y but it might be both? MOST IMPORTANT THING TO FIGURE OUT!
###Regression function needs review 

## Create a compiler file to run all these functions 



function initialize_coef(R:Results) #See PS5 for choices 
    @unpack a0, b0, a1, b1 = R
    a0 = 0.095
    b0 = 0.085
    b1 = 0.999
    a1 = 0.999
    return a0, b0, a1, b1 
end

function SimulateCapitalPath(R::Results, P::Params, Epsilon,Z)
    @unpack pf_k = R
    @unpack N, T, burn = P

    kappa::Array{Float64,2} = zeros(T-1,2) #since the first one is Kbar2
    V::Array{Float64,2} = zeros{N,T} #Storing the panel. I hate that we have to do this column by column but rows come first 
    
    #do initial column seperately
    for n in 1:N
        V[n,1] = pf_k[11.55,Epsilon[n,1],11.55,Z[1]] #11.55 is a hardfix, should replace with steady state later 
    end
    kappa[1,1] = sum(V[:,1])/N
    kappa[1,2] = Z[1]
    
    for t in 2:T
        for n in 1:N
            V[n,t] = pf_k[V[n,t-1],Epsilon[n,t],kappa[t-1],Z[t]] #need to flip my epsilon I think 
        end
        kappa[t,1] = sum(V[t,:])/N
        kappa[t,2] = Z[t] #for the sort lateron
    end

    ##I might want a third column in kappa for either the previous or next average K, to act as y since I need that before I sort. 

    #burn the first 1000 columns/rows (kappa is weird)
    kappa_b::Array{Float64,2} = kappa[burn+1:T,:]
    V_b::Array{Float64,2} = V[:,burn+1:T]

    return kappa_b, V_b  
end

function state_sort(array) #figure this might be worth doing, will need to check dimensions
    n_b = length(array[:,1])
    good::Array{Float64,1} = [0]
    bad::Array{Float64,1} = [0]

    for t in 1:n_b #there must be a less janky way to do this but Julia does not like conditional deletes
        if array[t,2] = 1
            good = append!(good,array[t,1])
        else
            bad = append!(bad,array[t,1])
        end
    end
    l_g = length(good)
    l_b = length(bad)
    good = good[2:l_g] #remove the first line
    bad = bad[2:l_b]

   ##Might need to return a Nx2 matrix instead, in which case append won't work, will need to make a state_y array instead
    return good, bad 
end 

function log_autoreg(x,y,cons::Bool = true) #this one needs checking, especially how I get y 
    if cons == true
        n_x = length(x)
        lx = log.(x)
        xi = hcat(ones(n_x),lx) #create a new matrix with a column of ones, removes the last row because auto regressive
        xit = transpose(xi)
        coefs = (xit*xi)^(-1)*(xit*y) #linear regression 
        intercept = coefs[1]
        coef = coefs[2]
        error = y - xi*coefs
        y_m = y-ones(n_x)*(ones(n_x)'*ones(n_x))^(-1)*(ones(n_x)'*y)
        R2 = (transpose(error)*error)/(transpose(y_m)*y_m)
    return intercept, coef, R2
    else
        lx = log.(x) 
        xit = transpose(xi)
        coef = (xit*xi)^(-1)*(xit*y)
        error = y - xi*coef
        y_m = y-ones(n_x)*(ones(n_x)'*ones(n_x))^(-1)*(ones(n_x)'*y)
        R2 = (transpose(error)*error)/(transpose(y_m)*y_m)
    return coef, R2
    end
end


function EstimateRegression(kappa,epsilon,cons::Bool = true )
    #first do a sort based on z
    kappa_a, kappa_b = state_sort(kappa) #a is for good state

    #Might need to create a Y here or I can make it in the log_autoreg fn

    #Next run regressions
    ahat0, ahat1, R2[1] = log_autoreg(kappa_a,epsilon,cons) #not sure what the right Y is, it's not simply epsilon, it might be Kappa again
    ahat0, ahat1, R2[2] = log_autoreg(kappa_b,epsilon,cons) #not sure what the right Y is

    return ahat0, ahat1, bhat0, bhat1, R2
end 




function Update_Coef(R::Results,lambda)
    @unpack a0,b0,a1,b1, ahat0,ahat1,bhat0,bhat1 = R
    a0 = lambda .*ahat0 .+ (1 .-lambda) .*a0
    b0 = lambda .*bhat0 .+ (1 .-lambda) .*b0
    a0 = lambda .*ahat1 .+ (1 .-lambda) .*a1
    a0 = lambda .*bhat1 .+ (1 .-lambda) .*b1
    return a0,b0,a1,b1
end

function Solve_KS(P::Params, G::Grids, S::Shocks, R::Results, lambda::Float64=0.5, tol_up = 1e-3, m_gof::Float64 =1-1e-2, kill = 10000)
   #First step is to draw shocks
    Epsilon, Z = DrawShocks(S,P.N,P.T) #use the formula from the helpful functions to get our shocks
    #Second step is to set initial a0,b0,a1,b1

    R.a0, R.b0, R.a1, R.b1 = initialize_coef(R)

    stop = 0
    count = 0
    while stop = 0 && minimum(R.R2) < m_gof ##wrong way to do the goodness of fit loop 
        count +=1
        R.pf_k, R.pf_V = Bellman(P,G,S,R) #note this returns K' first, different than how we normally do it

        kappa, V_b = SimulateCapitalPath(R,P,Epsilon,Z) #I think this is all correct now 
        R.ahat0, R.ahat1, R.bhat0, R.bhat1, R.R2 = EstimateRegression(kappa) #Still need to write this 

        if count < kill && (abs(R.ahat0-R.a0)+abs(R.ahat1-R.a1)+abs(R.bhat0-R.b0)+abs(R.bhat1-R.b1)) > tol_up
            R.a0, R.b0, R.a1, R.b1 = Update_Coef(R::Results, lambda) #Still need to write this up
        else 
            stop = 1
            println("The coefficients have converged")
        end
    end
end



