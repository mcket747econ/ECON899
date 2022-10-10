#So the Helpful functions file has drawshocks, a bellman with interpolation and our parameter structutres
#Need to remember to bring in both it and this model file in our compute code
#Beyond the usual packages, need "Optim" and "Interpolations" otherwise the Bellman won't work 


#To Do:  
##Full model solver
###Establish structutre
###Put in precomposed functions
###Figure out how to do the goodness of fit. Seems that we might want an outer while loop. 

##Write SimulateCapitalPath
###Ultimate goals, get average capital of each z, return a 1xT-1 vector
###Do one column at a time
###Burn Kappa before returning 

## write EstimateRegression
###Make sure there are two equations, one for each state of the world
###Will need to sort data across states
###Returns 4 estimated coefficients 



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

    kappa::Array{Float64,1} = zeros(T-1) #since the first one is Kbar2
    V::Array{Float64,2} = zeros{T,N} #Storing. Julia indexes columns first, right? 
    
    #do initial column seperately
    for n in 1:N
        V[1,n] = pf_k[11.55,Epsilon[t,n],11.55,Z[t]] #not sure what z_t in the notes refers to, 11.55 is a hardfix, should replace with steady state later 
    end
    kappa[1] = sum(V[1,:])/N
    
    for t in 2:T
        for n in 1:N
            V[t,n] = pf_k[V[t-1,n],Epsilon[t,n],kappa[t-1],Z[t]] #not sure what z_t in the notes refers to 
        end
        kappa[t] = sum(V[t,:])/N
    end

    #burn the first 1000 columns 
    kappa_b::Array{Float64,1} = kappa[burn+1:T]
    V_b::Array{Float64,2} = V[burn+1:T,:]

    return kappa_b, V_b  
end


function EstimateRegression()
    #first do a sort based on z

    #Next run regressions, don't forget column vector of ones (or find a package to do this)


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
    Z , Epsilon = DrawShocks(S,P.N,P.T) #use the formula from the helpful functions to get our shocks
    #Second step is to set initial a0,b0,a1,b1

    R.a0, R.b0, R.a1, R.b1 = initialize_coef(R)

    stop = 0
    count = 0
    while stop = 0 && R.R2 < m_gof ##wrong way to do the goodness of fit loop 
        count +=1
        R.pf_k, R.pf_V = Bellman(P,G,S,R) #note this returns K' first, different than how we normally do it

        kappa, V_b = SimulateCapitalPath(R,P,Epsilon,Z) #Still need to write this 
        R.ahat0, R.ahat1, R.bhat0, R.bhat1, R.R2 = EstimateRegression(kappa) #Still need to write this 

        if count < kill && (abs(R.ahat0-R.a0)+abs(R.ahat1-R.a1)+abs(R.bhat0-R.b0)+abs(R.bhat1-R.b1)) > tol_up
            R.a0, R.b0, R.a1, R.b1 = Update_Coef(R::Results, lambda) #Still need to write this up
        else 
            stop = 1
            println("The coefficients have converged")
        end
    end
end



