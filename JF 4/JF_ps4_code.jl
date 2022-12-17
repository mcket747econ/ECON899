using DataFrames, Random, Parameters, Distributions, Accessors, CSV
using Pkg;Pkg.add("CSV")
@with_kw struct params
    alpha::Float64 = 2
    lambda::Int64 = -4
    beta::Float64 = 0.99
    i_bar::Int64 = 8
    ps::Int64 = 1
    pr::Int64 = 4
    c_0::Int64 = 0
    c_1::Int64 = 1
    c_prob::Float64 = 0.5
    



end


S = DataFrame(CSV.File("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS4/PS4/PS4_state_space.csv"))
F_a0 = DataFrame(CSV.File("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS4/PS4/PS4_transition_a0.csv"))
F_a1 = DataFrame(CSV.File("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS4/PS4/PS4_transition_a1.csv"))
sim_data = DataFrame(CSV.File("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS4/PS4/PS4_simdata.csv"))

@with_kw mutable struct results
    exp_val_func::Array{Float64,3}
    val_func_0::Array{Float64,3}
    val_func_1::Array{Float64,3}










end 

function initialize(S)
    P = params()
    exp_val_func = zeros(size(S)[1],size(S)[1],size(S)[1])
    val_func_0 = zeros(size(S)[1],size(S)[1],size(S)[1])
    val_func_1 = zeros(size(S)[1],size(S)[1],size(S)[1])





    R = results(exp_val_func,val_func_0,val_func_1)
    return P, R
end
P,R = initialize(S)

function payoff(P::params, R::results,a,i,c,p)
    if a == 0 && i ==0
        payoff = P.lambda_1*c 
    elseif a == 0 && i > 0
        payoff = P.alpha*c 
    else
        payoff = P.alpha*c - p 
    end
    return payoff 
end


function value_func(P::params, R::results,i,c,p,S)
    val_func = zeros(size(S)[1],size(S)[1],size(S)[1])
    val_0 = zeros(size(S)[1],size(S)[1],size(S)[1])
    val_1 = zeros(size(S)[1],size(S)[1],size(S)[1])                                                                                         
    for ind = 1:size(S)[1]
        for a = 0:1
            v = P.alpha*S[ind,"C"] - S[ind,"P"]
            println(S[ind,"I"] + a - S[ind,"C"])
            #i_p = minimum(P.i_bar,(S[ind,"I"] + a - S[ind,"C"]))
            if S[ind,"I"] + a - S[ind,"C"] > P.i_bar
                i_p = i_bar 
            else
                i_p = S[ind,"I"] + a - S[ind,"C"]
            end 

            if a == 0
            val_0 = payoff(P,R,a,S[ind,"I"],S[ind,"C"],S[ind,"P"]) + P.beta*0.5*(val_func[i_p,0,P.pr]*0.9 + val_func[i_p,0,P.ps]*0.1)
                                                                                                     + P.beta*0.5*(val_func[i_p,1,P.pr]*0.9 + val_func[i_p,1,P.ps]*0.1)
            elseif a == 1
                val_1 = payoff(P,R,a,S[ind,"I"],S[ind,"C"],S[ind,"P"]) + P.beta*0.5*(val_func[i_p,0,P.pr]*0.9 + val_func[i_p,0,P.ps]*0.1)
                + P.beta*0.5*(val_func[i_p,1,P.pr]*0.9 + val_func[i_p,1,P.ps]*0.1)
            end
            if val_0 > val_ 1
                val = val_0
            else
                val = val_1
            end               
            return val,val_1,val_0
        end
        val_func[S[ind,"I"],S[ind,"C"],S[ind,"P"]] = val
    end

     return val_func,val_0,val_1                                                                                           
end

R.val_func,R.val_func_0, R.val_func_1 = value_func(P,R,S[:,"I"],S[:,"C"],S[:,"P"],S)

function expected_val_func(P,R,i,p,c,s)
    val_func_0 = R.val_func_0
    val_func_1 = R.val_func_1
    v_tilde = val_func_1 - val_func_0
    p_1s = exp.(v_tilde./(1 .+ exp.(v_tilde)))
    p_0s = 1 .- p_1s 
    e_0 = P.gamma .- log.(p_0s)
    e_1 = P.gamma .- log.(p_1s)
    exp_val_func = p_0s*(val_func_0 .+ e_0) + p_1s*(val_fun_1 .+ e_1)



    return exp_val_func
end



