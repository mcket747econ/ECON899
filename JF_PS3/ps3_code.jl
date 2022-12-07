# using Parameters, DataFrames, StatFiles, Random, Distributions, Optim, LinearAlgebra

cd("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3")

ind_characteristics = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Simulated_type_distribution.dta"))
car_ch = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Car_demand_characteristics_spec1.dta"))
iv_specification = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Car_demand_iv_spec1.dta"))






@with_kw struct Param
    T::Int64 = length(unique(car_ch[:,"Year"]))  ##The number of markets(Years)
    # delta_iia::Array{Float64,1} = Vector(car_ch[:,"delta_iia"]) ##the iia delta values
    simdist::Array{Float64,1} = Vector(ind_characteristics[:,"Var1"]) ##Vector containing the simulated individual characteristics
    mz::Array{Float64,1} = Vector(car_ch[:,"price"]) ##Vector containing prices for all models
    lambda::Float64 = 0.6
    eps1::Float64 = 1 
    eps0::Float64 = 10e-12
    years::Array{Float64,1}= sort(unique(car_ch[:,"Year"]))
    share_ch::Array{Float64,1} = Vector(car_ch[:,"share"])








end

@with_kw mutable struct Results  ##Include elements that will change over the course of the iterations
    produc_t = Array{Array{Float64,1},1}()
    delta_iia::Array{Float64,1} 
    #res_norm::Array{Float64,1} 






end 


function Initialize() #Initialize parameter struct and results struct
    P = Param()
    produc_t = Array{Array{Float64,1},1}()
    delta_iia = Vector(car_ch[:,"delta_iia"]) ##the iia delta values
    R = Results(produc_t,delta_iia)

    return P ,R
end 

P,R = Initialize()



function precompute(carch_ch_mt,years,produc_t)  ##Precompute the productids for each particular year
    produc_1 = 0
    for i = 1:31
        produc_1 = findall(x->x==years[i],carch_ch_mt[:,2])
        push!(produc_t,produc_1)

    end
    return produc_t
end


function individual_char(P::Param, produc_t)  
    Zf = Array{Array{Float64,2},1}() ##Structure that will contain all the prices associated with products in a particular market(Year)
    Z = zeros(P.T,1)
    for i=1:P.T
        for j=1:1
            Z = P.mz[Int64.(produc_t[i])]*P.simdist' 
            push!(Zf,Z)
            #Zf[i,j] = Z

        end
    end
    return Zf
end



function value(lambda,t,Z)   ##Evaluating the value of individual characteristics,and prices u
    mu = Array{Float64,2}
    l = 1 ##Number of nonlinear attributes
    mu = lambda.*Z[t]
    mu_exp = exp.(mu)
    return mu_exp
end

function demand(mu,t,delta,produc_t) ##Creating the Demand Function
    ev = exp.(delta[Int64.(produc_t[t])]).*mu
    #ev = exp.(delta[Int64.(produc_t[t])] .- logsumexp.(delta[Int64.(produc_t[t])])).*mu
    market_share = ev./(1 .+ ev)
    #println("market share is ", market_share)
    shat = mean(market_share,dims=2)
    return market_share, shat

end


function jacobian(ms)

    square = ms*(1 .- ms)'
    mat2 = ms*ms'#./100
    sz = size(mat2)
    jacob = mean.(eachcol(square)) .- (1/100).*mat2 + Diagonal(zeros(size(mat2)))  #(1/100).*(1 - I).*mat2
    return jacob
end





function inverse(P::Param,R::Results,eps0,eps1,share_ch,T_thisrun,shat,delta_iia,mu,produc_t,Z)
    for t = 1:T_thisrun
        value(P.lambda,t,Z) ##Calci
        rowid = produc_t[t]
        f = 100
        res_norm::Array{Float64,1} = [0]
        i = 1
        jacob = Matrix{Float64}
        while norm(f) > eps1   ##While loop for when norm is greater than 1
            #jacob = 0
            #println("norm f:", norm(f))
            ms, shat= demand(mu,T_thisrun,delta_iia,produc_t)
            #println(i)
            #println("shat is",shat)
            jacob = jacobian(ms)
            # println(P.share_ch[Int64.(rowid)])
            # println(log.(shat))
            # println(log.(P.share_ch[Int64.(rowid)]))
            f = log.(P.share_ch[Int64.(rowid)]) .- log.(shat)
            # print(f)
            R.delta_iia[Int64.(rowid)] = R.delta_iia[Int64.(rowid)] .+ f
            println("This is R Delta", R.delta_iia[Int64.(rowid)])
            #println(norm(f))
            append!(res_norm,norm(f))
            #print("res norm really is ", res_norm)
            i += 1
        end 

        while norm(f) > eps0  ##While loop for when norm is less than 1 but above 10e-12
            #jacobian = 1
            println("norm f:", norm(f))
           #println("this is delt ", delta_iia[Int64.(produc_t[T_thisrun])])  
            ev = (delta_iia[Int64.(produc_t[T_thisrun])]).*mu  
           # println("this is delt ", delta_iia[Int64.(produc_t[T_thisrun])])  
            # println("this is ev ", ev)   
            # println("share is", (ev)./(1 .+ ev))
            # println("market share is", mean((ev)./(1 .+ ev),dims=2))
            ms,shat = demand(mu,T_thisrun,delta_iia,produc_t)
            #  ms = (ev)./(1 .+ ev)
            #  shat = mean((ev)./(1 .+ ev),dims=2)
            # println("market share2 is", ms)
           #  println("this is ms ", ms)   
            # println("shat2 is ", shat)
            f = log.(P.share_ch[Int64.(rowid)]).-log.(shat)
            jacob = jacobian(ms)
            # println(" Why is this Na ", R.delta_iia[Int64.(rowid)])
            # println(" Why is this f Na ", f)
            #println(jacob./shat)
            R.delta_iia[Int64.(rowid)] = R.delta_iia[Int64.(rowid)]  + inv(jacob./shat)*f
            println("This is R Delta", R.delta_iia[Int64.(rowid)])
            #println(norm(f))
            append!(res_norm,norm(f))  
            #println("res norm really really is ", res_norm)    
        end     
        return res_norm   
    end
    #println("res norm is ", res_norm)
    return res_norm
end 










function routine_overall(P::Param,R::Results,t)  ##Overall code function
    P, R = Initialize()
    res_norm::Array{Float64,1} = [0]
    produc_t = Array{Array{Float64,1},1}()
    jacob = Matrix{Float64}
    Zf = Array{Array{Float64,2},1}()
    carch_ch_mt = Matrix(car_ch)
    produc_t = precompute(carch_ch_mt,P.years,produc_t)
    Zf = individual_char(P,produc_t)
    mu =value(P.lambda,t,Zf)
    market_share,shat  = demand(mu,t,R.delta_iia,produc_t)
    #println("shat is ", shat)
    jacob = jacobian(market_share)
    #print(jacob)
    res_norm_init = inverse(P,R,P.eps0,P.eps1, P.share_ch,t,shat,R.delta_iia,mu,produc_t,Zf)
   #print(res_norm_init)
   #println("res norm init ", res_norm_init)
    return res_norm_init  ##return vector of norms
end




# res_norm_init = inverse(P,R,P.eps0,P.eps1, P.share_ch,t,shat,R.delta_iia,mu,produc_t,Zf)










