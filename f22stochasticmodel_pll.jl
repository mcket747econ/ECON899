#keyword-enabled structure to hold model primitives
@everywhere using Pkg;Pkg.add(["Plots","Parameters"])
@everywhere using Parameters, Plots, SharedArrays
@everywhere @with_kw struct Primitives

    β::Float64 = 0.99 #discount rate
    δ::Float64 = 0.025 #depreciation rate
    α::Float64 = 0.36 #capital share
    k_min::Float64 = 1 #capital lower bound
    k_max::Float64 = 60.0 #capital upper bound
    k_grid::Array{Float64,1} = collect(range(1.0, length = 1000, stop = 60.0)) #capital grid
    markov::Array{Float64,2} = [0.977 0.023; 0.074 0.926] #Transition Matrix
    z_grid::Array{Float64,1} =[1.25,0.20] #Possible Values of z
    nk::Int64 = length(k_grid) #number of capital grid points
    nz::Int64 = length(z_grid)#The number of productivity states
end

#structure that holds model results

@everywhere mutable struct Results
    val_func::SharedArray{Float64,2} #value function
    pol_func::SharedArray{Float64,2} #policy function
end

#function for initializing model primitives and results
@everywhere function Initialize()
    prim = Primitives() #initialize primtiives
    val_func = zeros(prim.nk, prim.nz) #initial value function guess
    pol_func = zeros(prim.nk, prim.nz) #initial policy function guess
    res = Results(val_func, pol_func) #initialize results struct
    V_iterate(prim, res)
    prim, res #return deliverables
end
#iterate
@everywhere function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-3)
    error = 100
    n =  0 #counter
    while error>tol #begin iteration
        n+=1
        v_next = Bellman(prim, res) #spit out new vectors
        error = maximum(abs.(v_next-res.val_func)) #reset error level
        res.val_func = v_next #update value function
    end
    println("Value functions converged in ", n, " iterations.")
end

#Bellman Operator
@everywhere function Bellman(prim::Primitives,res::Results)
    @unpack k_grid, β, δ, α, nk, nz, z_grid, markov = prim #unpack model primitives
    v_next = zeros(nk, nz) #next guess of value function to fill

    #for exploiting monotonicity of policy function; can't do this in parallel
    @sync @distributed for i in nz*nk #single loop over entire state space, as discussed in TA section. Yes this takes heavily from the example
        i_z = mod(i,nz) + 1
        i_k = mod(ceil(Int64,i/nz),nk) + 1 
        k = k_grid[i_k] #value of k today
        z = z_grid[i_z] #getting the proper state

        Prob = markov[i_z,:]

        y = z*k^α  #need to define outside loop

        #I didn't want to take this from the example but otherwise am not sure how to best speed up updating and make the right choice
        v_0 = log(0)
        c_0 = log(0)
        k
    
        for (kp_index, k_temp) in enumerate(k_grid) #loop over possible selections of k', exploiting monotonicity of policy function
            c = y + (1-δ)*k - k_temp #consumption given k' selection
            v_prime = Prob[1]*res.val_func[kp_index,1] + Prob[2]*res.val_func[kp_index,2]
            if c<0 #check for positivity
                val = log(0) + β*v_prime#compute value
            else 
                val = log(c) + β*v_prime
            end
            if val > v_0
                v_0 = val
                c_0 = c
            end
        end
        #updates 
        v_next[i_k,i_z] = v_0 #update value function
    end
    v_next #return next guess of value function
end

####Above is currently not working. I think it isn't distributing the loops right between the iterator and the bellman equation

#Value function iteration


#solve the model
#prim, res = Initialize() #initialize primitive and results structs


@everywhere function Solve_model(prim::Primitives, res::Results)
    V_iterate(prim, res) #in this case, all we have to do is the value function iteration!
end
##############################################################################
#@elapsed prim, res = Initialize() #solve the model.
#Plots.plot(prim.k_grid, res.val_func) #plot value function
#Plots.plot(prim.k_grid, res.pol_func) #plot value function