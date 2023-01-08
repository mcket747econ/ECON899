#using Pkg;Pkg.add("CSV") Load csv package in order to be able to load csv files

using DataFrames, Random, Parameters, Distributions, Accessors, CSV, LinearAlgebra, Optim

cd("/Users/jacobbills/Desktop/Economics/Econ 899/JF 4/")

S = DataFrame(CSV.File("./PS4/PS4_state_space.csv"))
F_a00 = DataFrame(CSV.File("./PS4/PS4_transition_a0.csv"))
F_a0 = Matrix(F_a00[:,3:38])
F_a11 = DataFrame(CSV.File("./PS4/PS4_transition_a1.csv"))
F_a1 = Matrix(F_a11[:,3:38])

sim_data = DataFrame(CSV.File("./PS4/PS4_simdata.csv"))

@with_kw struct params  ##
    alpha::Float64 = 2
    #lambda::Int64 = -4 #this will be modified, needs to be in a mutable structure
    beta::Float64 = 0.99
    gamma::Float64 = 0.57721 #Euler's Constant
    i_bar::Int64 = 8
    ps::Int64 = 1 #Sale Price
    pr::Int64 = 4   #Regular Price
    c_0::Int64 = 0 #0 consumption shock
    c_1::Int64 = 1 #non-zero consumption shock
    c_prob::Float64 = 0.5 #Probability of a particular consumption shock
    pr_prob::Float64 = 0.9 #Probability of a regular price 
    a::Array{Float64,1} = [0, 1] #choices of different actions
    na::Int64 = length(a) #number of actions
    p::Array{Float64,1} = [pr,ps] #Vector of prices
    np::Int64 = length(p) #number of price options
    c::Array{Float64,1} = [c_0, c_1] #Array of consumption shocks
    nc::Int64 = length(c) #Length of consumption shock array
    ns::Int64 = size(S)[1] #number of potential states
    i::Array{Float64,1} = [0,1,2,3,4,5,6,7,8] #Investment options array
    ni::Int64 = length(i) #number of investment options

end


@with_kw mutable struct results
    exp_val_func::Array{Float64,3}
    exp_val_vector::Array{Float64,1}
    exp_val_func_ccp::Array{Float64,1}
    val_func_0::Array{Float64,3}
    val_func_1::Array{Float64,3}
    val_func::Array{Float64,3}
    val_func_ccp::Array{Float64,1}
    p_hat::Array{Float64,1}
    lambda::Float64 = -4.0 #so that it is mutable 
end



function initialize(S)
    P = params()
    exp_val_func = zeros(size(S)[1],size(S)[1],size(S)[1])
    exp_val_vector = zeros(P.ns)
    exp_val_func_ccp = zeros(size(S)[1])
    val_func_0 = zeros(P.ni,P.nc,P.np)
    val_func_1 = zeros(P.ni,P.nc,P.np)
    val_func = zeros(P.ni,P.nc,P.np)
    val_func_ccp = zeros(P.ns)
    p_hat = fill(0.5,P.ns)
    lambda::Float64 = -4.0


    R = results(exp_val_func,exp_val_vector,exp_val_func_ccp,val_func_0,val_func_1,val_func,val_func_ccp,p_hat,lambda)
    return P, R
end
P,R = initialize(S)

function payoff(P::params, R::results,a,i,c,p)  ##Per Period Payoff Function 
    if a == 0 && i ==0   ##Payoff if don't restock and no inventory
        if c > 0   ##Positive consumption shock
            payoff = R.lambda ## PAy the stockout penalty
        else
            payoff = 0 #Otherwise no penalty 
        end
    elseif a == 0 && i > 0 #Payoff if we don't restock and our inventory is non zero
        payoff = P.alpha*c   
    elseif a == 1                    #If we choose to restock        
        payoff = P.alpha*c - p  #Payoff when choosing to restock 
    end
    return payoff 
end


function value_func(P::params, R::results,S)
    val_func = zeros(P.ni,P.nc,P.np)
    val_0 = 0
    val_1 = 0
    val0 = zeros(P.ni,P.nc,P.np) ##Choice 0 Specific Value Function
    val1 = zeros(P.ni,P.nc,P.np) #Choice 1 specific value function
    #val =  zeros(size(S)[1],size(S)[1],size(S)[1])                                                                                        
    for ind = 1:size(S)[1] ##Iterate over the index of each state row
        for a = 0:1 ##For each potential asset choice
            #v = P.alpha*S[ind,"C"] - S[ind,"P"] 
       #     #println(S[ind,"I"] + a - S[ind,"C"])
            #i_p = minimum(P.i_bar,(S[ind,"I"] + a - S[ind,"C"])),,,,,,,,,,,,,,,,,
            # Conditions for the investment value next period i_p
            # if S[ind,"I"] + a - S[ind,"C"] > P.i_bar ##
            #     i_p = P.i_bar ##i_p = i_bar if the amount of investment this period plus the restock decision minus the consumption decision is greater than i_bar
            # else
            #     i_p = S[ind,"I"] + a - S[ind,"C"] 
            # end 

            i_p = minimum([S[ind,"I"] + a - S[ind,"C"],P.i_bar])

            if a == 0  ##If we choose not to restock
                val_0 = payoff(P,R,a,S[ind,"I"],S[ind,"C"],S[ind,"P"])  ##Choice specific value function equals 
                #.5 is the probability of a particular consumption shock
                #0.9 is the probability of the price being P_r
                + P.beta*P.c_prob*((val_func[findall(x->x==i_p,P.i),findall(x->x==0,P.c),findall(x->x==4,P.p)])*P.pr_prob+ (val_func[findall(x->x==i_p,P.i),findall(x->x==0,P.c),findall(x->x==1,P.p)])*(1-P.pr_prob)) 
                + P.beta*P.c_prob*((val_func[findall(x->x==i_p,P.i),findall(x->x==1,P.c),findall(x->x==4,P.p)])*P.pr_prob + (val_func[findall(x->x==i_p,P.i),findall(x->x==1,P.c),findall(x->x==1,P.p)])*(1-P.pr_prob))
                 #R.val_func[i_p,1/2,.9*P.P_r +.1*P.P_s] = payoff(P,R,a,i_p,S[ind,"C"],S[ind,"P"]) + P.beta*0.5*(val_func[i_p,0,P.pr]*0.9 + val_func[i_p,0,P.ps]*0.1) + P.beta*0.5*(val_func[i_p,1,P.pr]*0.9 + val_func[i_p,1,P.ps]*0.1)
                val0[findall(x->x==S[ind,"I"],P.i),findall(x->x==S[ind,"C"],P.c),findall(x->x==S[ind,"P"],P.p)] .= val_0 #Set choice 0 specific value function
                #at a level of investment, consumption and price from this particular loop equal to val_0
            else 
                val_1 = payoff(P,R,a,S[ind,"I"],S[ind,"C"],S[ind,"P"]) 
                + P.beta*0.5*((val_func[findall(x->x==i_p,P.i),findall(x->x==0,P.c),findall(x->x==4,P.p)])*0.9 + (val_func[findall(x->x==i_p,P.i),findall(x->x==0,P.c),findall(x->x==1,P.p)])*0.1) 
                + P.beta*0.5*((val_func[findall(x->x==i_p,P.i),findall(x->x==1,P.c),findall(x->x==4,P.p)])*0.9 + (val_func[findall(x->x==i_p,P.i),findall(x->x==1,P.c),findall(x->x==1,P.p)])*0.1)
                val1[findall(x->x==S[ind,"I"],P.i),findall(x->x==S[ind,"C"],P.c),findall(x->x==S[ind,"P"],P.p)] .= val_1

            end
    
        end 
        if val_0 > val_1
            val_func[findall(x->x==S[ind,"I"],P.i),findall(x->x==S[ind,"C"],P.c),findall(x->x==S[ind,"P"],P.p)] .= val_0
        else
            val_func[findall(x->x==S[ind,"I"],P.i),findall(x->x==S[ind,"C"],P.c),findall(x->x==S[ind,"P"],P.p)] .= val_1
        end               
            #return val,val_1,val_0
    end

     return val_func,val0,val1   #Return the value function and the choice specific value functions                                                                              
end


#R.val_func,R.val_func_0,R.val_func_1 = value_func(P,R,S) #Assign value functions to their arrays in the results struct

function overall_iteration(P::params,R::results,S,tol::Float64=1e-5) ##This conducts our iteration in order to converge to the correct value
    new_val_func = zeros(P.ni,P.nc,P.np)
    val0 = zeros(P.ni,P.nc,P.np)
    val1 = zeros(P.ni,P.nc,P.np)
    error = 100
    n = 1
    while error > tol   ##Set our threshold 
        n+=1 #Counter
        new_val_func,val0,val1 = value_func(P,R,S) #Run our value function 
       # #println("error ", maximum(abs.(new_val_func - R.val_func))) #Print the eror result to keep track of what is happening
        error = maximum(abs.(new_val_func - R.val_func)) #Calculate the error between our new value function and the only one stored in R.val_func
        R.val_func = new_val_func ## Add the new value function from the current iteration since the tolerance between it and the old one was bigger than our threshold. 
        R.val_func_0 = val0
        R.val_func_1 = val1
        #println("Iteration ", n) #Print the iteration value to keep track of what is happening
        
    end
    ##println("Value functions converged in ", n, " iterations.")
end




function expected_val_func(P::params,R::results)     ##Calculating the expected value function by weighting between the choice specific value functions
    exp_val_func = zeros(P.ni,P.nc,P.np)
    exp_val_vector = zeros(P.ns)

    val_func_0 = R.val_func_0
    val_func_1 = R.val_func_1
    v_tilde = R.val_func_1 .- val_func_0  ##Per the formula, subtract val0 from val1 to yield val tilde
    p_1s = (1 .+ exp.(-1 .* v_tilde)).^(-1)
    #p_1s = exp.(v_tilde./(1 .+ exp.(v_tilde)))  #Not sure why this version produces erroneous values. It may be something subtle. 
   # println(p_1s) # 
    p_0s = 1 .- p_1s 
    
    e_0 = P.gamma .- log.(p_0s)
    e_1 = P.gamma .- log.(p_1s)
    exp_val_func .= p_0s.*(val_func_0 .+ e_0) + p_1s.*(val_func_1 .+ e_1)
    #exp_val_func = log()
    exp_val_vector[1:9] = exp_val_func[:,1,1]
    exp_val_vector[10:18] = exp_val_func[:,2,1]
    exp_val_vector[19:27] = exp_val_func[:,1,2]
    exp_val_vector[28:36] = exp_val_func[:,2,2]



    return exp_val_func,exp_val_vector
end

function CCP(P::params,R::results)
    p_hat = zeros(P.ns)
    stat_id = unique(sim_data[:,"state_id"]) #Get a vector of the unique values of state_id column. Aka, get all unique state_ids. 
    length_stat = length(stat_id) ##Get the length of this vector so that we know the max value over which to iterate. 
    for stat in 1:length_stat ##For loop over unique states
        ##There is another way to do this, which would be looping over each row in the simulation matrix(5000 values)
        #But this way seemed more parsimonious
       numb_in_state = size(sim_data[sim_data.state_id.==stat_id[stat],:])[1]
       #get the number of rows where the the state is the curent valuer of stat_id
        #for row in length(sim_data[sim_data.state_id.==stat,"choice"])
        number_restocking = size(sim_data[(sim_data.state_id.==stat_id[stat]).&& (sim_data.choice.==1),:])[1]
        p_hat[stat] = number_restocking/numb_in_state
        println(p_hat[stat], stat)
        #numerator: the size of the subset of the simulation matrix where the state is the curent value of stat_id and 
        #the value of the choice column is 1(i.e., where the state is state_id[stat] and the individual chooses to restock)
        #denominator: the number of values that are in the state stat_id[stat]
        #phat: numerator/denominator. The share of individuals in state state_id[stat] who choose to restock
    end 
    ####Constrain p_hat values 
    for i = 1:length(p_hat)
        if p_hat[i] < 0.001
            p_hat[i] = 0.001
        elseif p_hat[i] > 0.999
            p_hat[i] = 0.999
        end 
    end 

    return p_hat
end 

function payoff_ccp(P::params, R::results,S)
    val_func = zeros(P.ni,P.nc,P.np)
    val_0 = 0
    val_1 = 0
    val0 = zeros(P.ns) ##Choice 0 Specific Value Function
    val1 = zeros(P.ns) #Choice 1 specific value function
    #val =  zeros(size(S)[1],size(S)[1],size(S)[1])                                                                                        
    for ind = 1:size(S)[1] ##Iterate over the index of each state row
        for a = 0:1 ##For each potential asset choice
            v = P.alpha*S[ind,"C"] - S[ind,"P"] 
          #  println(S[ind,"I"] + a - S[ind,"C"])
           # println(S[ind,"I"] + a - S[ind,"C"])
            #i_p = minimum(P.i_bar,(S[ind,"I"] + a - S[ind,"C"])),,,,,,,,,,,,,,,,,
            # Conditions for the investment value next period i_p
            if S[ind,"I"] + a - S[ind,"C"] > P.i_bar ##
                i_p = P.i_bar ##i_p = i_bar if the amount of investment this period plus the restock decision minus the consumption decision is greater than i_bar
            else
                i_p = S[ind,"I"] + a - S[ind,"C"] 
            end 

            if a == 0  ##If we choose not to restock
                val_0 = payoff(P,R,a,S[ind,"I"],S[ind,"C"],S[ind,"P"])  ##Choice specific value function equals 
                #.5 is the probability of a particular consumption shock
                #0.9 is the probability of the price being P_r
                #+ P.beta*P.c_prob*((val_func[findall(x->x==i_p,P.i),findall(x->x==0,P.c),findall(x->x==4,P.p)])*P.pr_prob+ (val_func[findall(x->x==i_p,P.i),findall(x->x==0,P.c),findall(x->x==1,P.p)])*(1-P.pr_prob)) 
                #+ P.beta*P.c_prob*((val_func[findall(x->x==i_p,P.i),findall(x->x==1,P.c),findall(x->x==4,P.p)])*P.pr_prob + (val_func[findall(x->x==i_p,P.i),findall(x->x==1,P.c),findall(x->x==1,P.p)])*(1-P.pr_prob))
                 #R.val_func[i_p,1/2,.9*P.P_r +.1*P.P_s] = payoff(P,R,a,i_p,S[ind,"C"],S[ind,"P"]) + P.beta*0.5*(val_func[i_p,0,P.pr]*0.9 + val_func[i_p,0,P.ps]*0.1) + P.beta*0.5*(val_func[i_p,1,P.pr]*0.9 + val_func[i_p,1,P.ps]*0.1)
                val0[ind] = val_0 #Set choice 0 specific value function
                #at a level of investment, consumption and price from this particular loop equal to val_0
            else a == 1
                val_1 = payoff(P,R,a,S[ind,"I"],S[ind,"C"],S[ind,"P"]) 
                #+ P.beta*0.5*((val_func[findall(x->x==i_p,P.i),findall(x->x==0,P.c),findall(x->x==4,P.p)])*0.9 + (val_func[findall(x->x==i_p,P.i),findall(x->x==0,P.c),findall(x->x==1,P.p)])*0.1) +
                #+ P.beta*0.5*((val_func[findall(x->x==i_p,P.i),findall(x->x==1,P.c),findall(x->x==4,P.p)])*0.9 + (val_func[findall(x->x==i_p,P.i),findall(x->x==1,P.c),findall(x->x==1,P.p)])*0.1)
                val1[ind] = val_1

            end
    
        end 
        if val_0 > val_1
            val_func[findall(x->x==S[ind,"I"],P.i),findall(x->x==S[ind,"C"],P.c),findall(x->x==S[ind,"P"],P.p)] .= val_0
        else
            val_func[findall(x->x==S[ind,"I"],P.i),findall(x->x==S[ind,"C"],P.c),findall(x->x==S[ind,"P"],P.p)] .= val_1
        end               
            #return val,val_1,val_0
    end

     return val_func,val0,val1   #Return the value function and the choice specific value functions                                                                              
end

function vbar_ccp(P::params,R::results,S,p_hat)
    #p_hat = CCP(P,R)
    payoff, payoff_0, payoff_1 = payoff_ccp(P,R,S)
    F = (F_a0).*(1 .- p_hat) .+ (F_a1).*(p_hat)
    F = Matrix(F)
   # println(F)
    e_0 = P.gamma .- log.(1 .- p_hat)
    e_1 = P.gamma .- log.(p_hat)
    ##println("size e",size(e_0))
    ##println("payoff_0",size(payoff_0))
    ##println("payoff",size(payoff))
    sz = size(P.beta*F)
    I_mat = Diagonal(ones(sz[1],sz[1]))
    inverse = inv(Matrix(I_mat .- P.beta*F))
    ##println("inverse", inverse)
    #inverse = Matrix(inverse)
  #  #println(length(p_hat))
   # #println(size(I_mat)[1])
   # #println(size("Size payoff and e", (payoff_0 .+ e_0)))
    weighted_v = ((1 .- p_hat).*(payoff_0 .+ e_0)) .+ p_hat.*(payoff_1 .+ e_1)
    #println(size(inverse))
    #println(size(weighted_v))
    V_barp = inverse*weighted_v
    #println("var_p", V_barp)





    return V_barp,payoff_0,payoff_1,F
end 


function overall_ccp_iteration(P::params, R::results, tol::Float64 = 10^(-10))
    exp_val_func_ccp = zeros(P.ni,P.nc,P.np)
    v_bar_new = zeros(P.ni,P.nc,P.np)
    payoff_0 = zeros(P.ns)
    payoff_1 = zeros(P.ns)
    new_phat = zeros(P.ns)
    F = zeros(P.ns,P.ns)
    v_tilde_ccp = zeros(P.ns)
    #payoff, payoff_0, payoff_1 = payoff_ccp(P,R)
    error = 100
    error2 = 100
    n = 1
    coef = fill(1/36,36)
    while error > 10^(-2)   ##Set our threshold 
        n+=1 #Counter
       
        exp_val_func_ccp,payoff_0, payoff_1, F = vbar_ccp(P,R,S,R.p_hat) #Run our value function 
        error =norm(exp_val_func_ccp .- R.exp_val_func_ccp)
        R.exp_val_func_ccp = exp_val_func_ccp
        # mf = F_a0.*R.p_hat+ F_a1.*(R.p_hat)
        # new_coef = mf'coef
        # error = norm(new_coef .- coef)
        # coef = new_coef
        # println(size(mf))    
        # #  R.exp_val_func_ccp .= v_bar_new
        # for i = 1:36
        #     R.exp_val_func_ccp[i] = mf[i,i]
        # end
        
    end

    while error2 > tol .&& n < 100
        n+=1
        ##println("error ", maximum(abs.(new_val_func_ccp - R.val_func_ccp))) #Print the eror result to keep track of what is happening
        #println("F_a1",F_a1)

        #v_tilde_ccp = (P.beta*(F_a1*exp_val_func_ccp) .+ ((payoff_1))) .- ((payoff_0) .+ P.beta*(F_a0*exp_val_func_ccp))
        #println("v_tilde_ccp ",v_tilde_ccp)
        #new_phat = (1 .+ exp.((-1)*v_tilde_ccp))
        exp_val_func_ccp,payoff_0, payoff_1, F = vbar_ccp(P,R,S,R.p_hat) #Run our value function 

        new_phat = (P.beta*(F_a1*exp_val_func_ccp) .+ ((payoff_1)))./((P.beta*(F_a1*exp_val_func_ccp) .+ ((payoff_1))) .+ ((payoff_0) .+ P.beta*(F_a0*exp_val_func_ccp)))
        new_phat = [if new_phat[x] > 0 && new_phat[x] < 1 new_phat[x] elseif new_phat[x] >= 1 .999 else .001 end for x=1:36] #constraining the frequency
        # mat_test = [payoff_0 payoff_1]
        # exp_test = exp.(mat_test)
        # sumr = sum(exp_test,dims=2)
        # new_phat = exp.(mat_test[:,2])./sumr

        #println(new_phat)
        #new_phat = Array(new_phat)
     #   #println(new_phat)
        ##println('size_phat', size(new_phat))
        ##println("error ", maximum(abs.(new_phat - R.p_hat))) #Print the eror result to keep track of what is happening
      #  #println(new_phat)
       # error = maximum(abs.(new_phat .- R.p_hat))
        error2 = norm( R.p_hat .- new_phat)
        #error2 = norm(new_phat .- R.p_hat) 
        #println(error2)
        #println(size(new_phat))
        ##println(new_phat)
        R.p_hat .= new_phat
        #R.exp_val_func_ccp = exp_val_func_ccp

        #println("Iteration ", n) #Print the iteration value to keep track of what is happening
        
    end
   println("Value functions converged in ", n, " iterations.")
   return R.p_hat
end

function inner(P::params, R::results,tol1::Float64= 10e-2,tol2::Float64=10e-10)
    overall_iteration(P,R,S,tol1)
    choice_prob = overall_ccp_iteration(P,R,tol2)
    return choice_prob
end


function ll(P::params, R::results,l,tol1::Float64= 10e-2,tol2::Float64=10e-10) #generates the (partial) log-likliehood 
    R.lambda = l[1]
    e_prob = inner(P,R,tol1,tol2)
    println("Done computing probabilities: lambda= ",l)
    partial_ll = 0 
    for i=1:5000
        x = sim_data[i,"state_id"] +1
        y = sim_data[i,"choice"]
        rel_p = [ if e_prob[x] >0 e_prob[x] else 10e-31 end][1]
        partial_ll = partial_ll + log(rel_p)*y + log(1-rel_p)*(1-y)
    end
    return -partial_ll 
end #this works since the partial log-likelihood (the sum of the log probability of each state) is consistent




function nfxp(P::params, R::results,data,l,tol1::Float64= 10e-2,tol2::Float64=10e-10) #finds the maximum likelihood 
    P,R = initialize(data) #re-initialize the data I guess
    #inner formula-> used to get the partial log-likelihood at each possible lambda 
    obj(lam) = ll(P,R,lam,tol1,tol2)
    #outer forumla-> does a search over possible lambdas 
    opt = optimize(obj, [l], BFGS())
    res = Optim.minimizer(opt)
    return res 
end

