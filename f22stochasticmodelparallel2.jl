#Original Code written in 2021, Modified by Matt McKetty Fall 2022
@everywhere using SharedArrays #must use shared arrays in order to allow for arrays to be updated while using parallel processing
@everywhere using Pkg;Pkg.add(["Plots","Parameters"])
@everywhere using Parameters, Plots #Load Necessary Packages
@everywhere @with_kw struct Primitives #Intialize model primitives

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

    #for exploiting monotonicity of policy function
    @sync @distributed  for k_index = 1:nk*nz #Use @sync and @distributed macro to ensure processing is spread out among multiple processors on computer
        candidate_max = -1e10#bad candidate max
        z_index = mod(k_index,nz) + 1 #index over good and bad states
        k_index =  mod(ceil(Int64, k_index/nz), nk) + 1 #index over values of k
        k = k_grid[k_index] #value of k
        z= z_grid[z_index] #values of z
        budget = z*k^α + (1-δ)*k #budget

        for kp_index in 1:nk #loop over possible selections of k', exploiting monotonicity of policy function
            c = budget - k_grid[kp_index] #consumption given k' selection
            if c>0 #check for positivity
                val = log(c) + β*sum(res.val_func[kp_index,:].*markov[z_index, :])#compute new value, sum of log and dot product of possible vs and their probabilities
                if val>candidate_max #check for new max value
                    candidate_max = val #update max value
                    res.pol_func[k_index,z_index] = k_grid[kp_index] #update policy function #update lowest possible choice
                end
            end
        end
        v_next[k_index,z_index] = candidate_max #update value function
    end
    v_next #return next guess of value function
end

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
