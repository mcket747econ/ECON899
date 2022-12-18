#using Pkg;Pkg.add("CSV") Load csv package in order to be able to load csv files

using DataFrames, Random, Parameters, Distributions, Accessors, CSV

S = DataFrame(CSV.File("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS4/PS4/PS4_state_space.csv"))
F_a0 = DataFrame(CSV.File("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS4/PS4/PS4_transition_a0.csv"))
F_a1 = DataFrame(CSV.File("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS4/PS4/PS4_transition_a1.csv"))
sim_data = DataFrame(CSV.File("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS4/PS4/PS4_simdata.csv"))

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
    pr_prob::Float64 = 0.9
    a::Array{Float64,1} = [0, 1]
    na::Int64 = length(a)
    p::Array{Float64,1} = [pr,ps]
    np::Int64 = length(p)
    c::Array{Float64,1} = [c_0, c_1]
    nc::Int64 = length(c)
    ns::Int64 = size(S)[1]
    i::Array{Float64,1} = [0,1,2,3,4,5,6,7,8]
    ni::Int64 = length(i)

    



end


@with_kw mutable struct results
    exp_val_func::Array{Float64,3}
    val_func_0::Array{Float64,3}
    val_func_1::Array{Float64,3}
    val_func::Array{Float64,3}










end 

function initialize(S)
    P = params()
    exp_val_func = zeros(size(S)[1],size(S)[1],size(S)[1])
    val_func_0 = zeros(P.ni,P.nc,P.np)
    val_func_1 = zeros(P.ni,P.nc,P.np)
    val_func = zeros(P.ni,P.nc,P.np)


    R = results(exp_val_func,val_func_0,val_func_1,val_func)
    return P, R
end
P,R = initialize(S)

function payoff(P::params, R::results,a,i,c,p)
    if a == 0 && i ==0
        if c > 0
            payoff = P.lambda*c 
        else
            payoff = 0
        end
    elseif a == 0 && i > 0
        payoff = P.alpha*c 
    else
        payoff = P.alpha*c - p 
    end
    return payoff 
end


function value_func(P::params, R::results,S)
    val_func = zeros(P.ni,P.nc,P.np)
    val_0 = 0
    val_1 = 0
    val0 = zeros(P.ni,P.nc,P.np)
    val1 = zeros(P.ni,P.nc,P.np)
    #val =  zeros(size(S)[1],size(S)[1],size(S)[1])                                                                                        
    for ind = 1:size(S)[1]
        for a = 0:1
            v = P.alpha*S[ind,"C"] - S[ind,"P"]
            println(S[ind,"I"] + a - S[ind,"C"])
            #i_p = minimum(P.i_bar,(S[ind,"I"] + a - S[ind,"C"])),,,,,,,,,,,,,,,,,
            if S[ind,"I"] + a - S[ind,"C"] > P.i_bar
                i_p = P.i_bar 
            else
                i_p = S[ind,"I"] + a - S[ind,"C"]
            end 

            if a == 0
                val_0 = payoff(P,R,a,S[ind,"I"],S[ind,"C"],S[ind,"P"]) 
                + P.beta*0.5*((val_func[findall(x->x==i_p,P.i),findall(x->x==0,P.c),findall(x->x==4,P.p)])*0.9 + (val_func[findall(x->x==i_p,P.i),findall(x->x==0,P.c),findall(x->x==1,P.p)])*0.1) 
                + P.beta*0.5*((val_func[findall(x->x==i_p,P.i),findall(x->x==1,P.c),findall(x->x==4,P.p)])*0.9 + (val_func[findall(x->x==i_p,P.i),findall(x->x==1,P.c),findall(x->x==1,P.p)])*0.1)
            #R.val_func[i_p,1/2,.9*P.P_r +.1*P.P_s] = payoff(P,R,a,i_p,S[ind,"C"],S[ind,"P"]) + P.beta*0.5*(val_func[i_p,0,P.pr]*0.9 + val_func[i_p,0,P.ps]*0.1) + P.beta*0.5*(val_func[i_p,1,P.pr]*0.9 + val_func[i_p,1,P.ps]*0.1)
                 val0[findall(x->x==S[ind,"I"],P.i),findall(x->x==S[ind,"C"],P.c),findall(x->x==S[ind,"P"],P.p)] .= val_0
            else a == 1
                val_1 = payoff(P,R,a,S[ind,"I"],S[ind,"C"],S[ind,"P"]) 
                + P.beta*0.5*((val_func[findall(x->x==i_p,P.i),findall(x->x==0,P.c),findall(x->x==4,P.p)])*0.9 + (val_func[findall(x->x==i_p,P.i),findall(x->x==0,P.c),findall(x->x==1,P.p)])*0.1) +
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

     return val_func,val0,val1                                                                                        
end

#R.val_func,R.val_func_0, R.val_func_1 = value_func(P,R,S[:,"I"],S[:,"C"],S[:,"P"],S)
R.val_func,R.val_func_0,R.val_func_1 = value_func(P,R,S)

function overall_iteration(P::params,R::results,S,tol::Float64=1e-5)
    new_val_func = zeros(P.ni,P.nc,P.np)
    val0 = zeros(P.ni,P.nc,P.np)
    val1 = zeros(P.ni,P.nc,P.np)
    error = 100
    n = 1
    while error > tol
        n+=1
        new_val_func,val0,val1 = value_func(P,R,S)
        println("error ", maximum(abs.(new_val_func - R.val_func)))
        println("error ", error)
        error = maximum(abs.(new_val_func - R.val_func))
        R.val_func = new_val_func
        R.val_func_0 = val0
        R.val_func_1 = val1
        println("Iteration ", n)
        
    end
    println("Value functions converged in ", n, " iterations.")
end

overall_iteration(P,R,S)



function expected_val_func(P::params,R::results)
    exp_val_func = zeros(P.ni,P.nc,P.np)
    val_func_0 = R.val_func_0
    val_func_1 = R.val_func_1
    v_tilde = R.val_func_1 .- val_func_0
    p_1s = (1 .+ exp.(-1 .* v_tilde)).^(-1)
    #p_1s = exp.(v_tilde./(1 .+ exp.(v_tilde)))
    println(p_1s)
    p_0s = 1 .- p_1s 
    
    e_0 = -4 .- log.(p_0s)
    e_1 = -4 .- log.(p_1s)
    exp_val_func .= p_0s.*(val_func_0 .+ e_0) + p_1s.*(val_func_1 .+ e_1)



    return exp_val_func
end

R.exp_val_func = expected_val_func(P,R)





