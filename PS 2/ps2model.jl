#keyword-enabled structure to hold model primitives
using Pkg;Pkg.add(["Plots","Parameters"])
using Parameters, Plots
@with_kw struct Primitives

    β::Float64 = 0.99.32 #discount rate
    α::Float64 = 1.5 #risk aversion
    a_min::Float64 = -2 #assets lower bound
    a_max::Float64 = 5 #assets upper bound
    a_grid::Array{Float64,1} = collect(range(-2, length = 1000, stop = 5)) #asset grid
    markov::Array{Float64,2} = [0.97 0.03; 0.5 0.5] #Transition Matrix
    s_grid::Array{Float64,1} =[1,0.5] #Possible Values of s
    na::Int64 = length(a_grid) #number of asset grid points
    ns::Int64 = length(s_grid)#The number of employment states
end

#structure that holds model results

mutable struct Results
    val_func::Array{Float64,2} #value function
    pol_func::Array{Float64,2} #policy function
end

mutable struct Distribution
    mu::Array{Float64,2} #distribution function
end

#function for initializing model primitives and results
function Initialize()
    prim = Primitives() #initialize primtiives
    val_func = zeros(prim.na, prim.ns) #initial value function guess
    pol_func = zeros(prim.na, prim.ns) #initial policy function guess
    res = Results(val_func, pol_func) #initialize results struct
    V_iterate(prim, res)
    prim, res #return deliverables
end

#q choices

q_lb = 0
q_ub = 1

#iterate
function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-3)
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
function Bellman(prim::Primitives,res::Results)
    @unpack a_grid, β, δ, α, na, ns, s_grid, markov = prim #unpack model primitives
    v_next = zeros(na, ns) #next guess of value function to fill

    #for exploiting monotonicity of policy function
    for a_index = 1:na, s_index=1:ns
        candidate_max = -1e10#bad candidate max
        a = a_grid[a_index] #value of k
        s= s_grid[s_index]
        budget = s+a-q*a_p  #budget

        for kp_index in 1:na #loop over possible selections of k', exploiting monotonicity of policy function
            c = budget - a_grid[kp_index] #consumption given k' selection
            if c>0 #check for positivity
                val = log(c) + β*sum(res.val_func[kp_index,:].*markov[s_index, :])#compute value
                if val>candidate_max #check for new max value
                    candidate_max = val #update max value
                    res.pol_func[a_index,s_index] = a_grid[kp_index] #update policy function #update lowest possible choice
                end
            end
        end
        v_next[a_index,s_index] = candidate_max #update value function
    end
    v_next #return next guess of value function
end

#Value function iteration


#solve the model
#prim, res = Initialize() #initialize primitive and results structs


function Solve_model(prim::Primitives, res::Results)
    V_iterate(prim, res) #in this case, all we have to do is the value function iteration!
end
##############################################################################
#@elapsed prim, res = Initialize() #solve the model.
#Plots.plot(prim.a_grid, res.val_func) #plot value function
#Plots.plot(prim.a_grid, res.pol_func) #plot value function

mu0 = ones(na,ns)
mu0 = mu0/na #initial value is uniformly distributed 


function t_star(prim::Primatives, res:Results)
    ##main idea: take the policy function and create Q, then iterate sum until done
    ##definitely not done yet 
    count = 1
    supnorm = 1
    while supnorm > .0001
        mu1 = zeros((na,ns))
        for si in 1:ns
            for ai in 1:na 
                Q[ai,si] =  
        end
    supnorm = norm((mu1-m0)/m0)
    count += 1
    mu0 = mu1
    end
end



function find_q(prim::Primatives, res::Results, dist::Distribution)
    mc_tol = 1
    q_ub = 1
    q_lb = 0
    while mc_tol > .0001
        q = (q_ub-q_lb)/2
        policy = Solve_model(prim,res,q)
        mu_star = t_star(policy, dist)
        ED = sum(policy*mu_star)    
        if ED > 0
            q_lb = q
        else if ED < 0
            q_ub = q
        end
        mc_tol = abs(mc_tol-ED)
    end
end


##gini index stuff

#need a function to create wealth distribution, the below isn't right. Want to know what the shape of mu is

function wealth(dist)
    for i in ns
     w(:,:,i)  = dist(:,:,i) + s_grid[i] #add earnings to assets
    end
    wealth_dist = vcat(w(:,:,1), w(:,:,2))
    return wealth_dist
end

#goal of the above is to turn the two state mu distribution into a single vector of assets plus earnings

function gini(dist)
    sorted = sort(dist)
    n = size(dist)
    coef = 2/n
    constant = (n+1)/n
    weighted_sum = sum([(i+1)*yi for i,yi in enumerate(sorted)])
    denom = sum(sorted)
    return coef*weighted_sum/(denom) - constant
end

##welfare functions

#need to figure out WFB 

function lambda(prim::Primitives,res::Results)
    @unpack a_grid, β, δ, α, na, ns, s_grid, markov = prim
    frac = 1/((1-α)*(1-β))
    for s in 1:ns
        for a in 1:na
            lam[a,s] = ((WFB + frac)/(val_func + frac))^(1/(1-α))-1
        end
    end
    return lam
end

function welfare(prim::Primitives, res::Results, lambda)
    WINC = sum(mu0*val_func)
    WG = sum(lambda*mu0)
    return [WINC,WG]
end
