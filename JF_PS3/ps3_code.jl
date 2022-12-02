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

@with_kw mutable struct Results
    produc_t = Array{Array{Float64,1},1}()
    delta_iia::Array{Float64,1} 
    #res_norm::Array{Float64,1} 






end 


function Initialize()
    P = Param()
    produc_t = Array{Array{Float64,1},1}()
    delta_iia = Vector(car_ch[:,"delta_iia"]) ##the iia delta values
    R = Results(produc_t,delta_iia)

    return P ,R
end 

P,R = Initialize()
# R.produc_t
# price_vector = zeros(6103,31 )
# z = zeros(31,1)
# share_ch = Vector(car_ch[:,"share"])

# carch_ch_mt = Matrix(car_ch)
# mz = Vector(car_ch[:,"price"])
# mz1 = Vector(car_ch[:,"price"])

# mz_test = Array{Float64,2}
# mz_test = zeros(6103,2)
# mz_test[:,1] = mz1
# mz_test[:,2] = x
# prices = car_p[xyz1]

Random.seed!(1234)
x = rand(Uniform(10,20),6103)
# xyz = getindex(carch_ch_mt[:,2],1997)
# xyz1 = findall(x->x==1985,carch_ch_mt[:,2])
# xyz2 = findall(x->x==1985,carch_ch_mt[:,2])
# xyz3 = findall(x->x==1985,carch_ch_mt[:,2])

# xy = [xyz1,xyz2,xyz3]
# xy[1]
# years = unique(car_ch[:,"Year"])

# Creating a product_id_code
function precompute(carch_ch_mt,years,produc_t)
    produc_1 = 0
    for i = 1:31
        produc_1 = findall(x->x==years[i],carch_ch_mt[:,2])
        push!(produc_t,produc_1)

    end
    return produc_t
end
# produc_t = Array{Array{Float64,1},1}()
# produc_t = precompute(carch_ch_mt,P.years,produc_t)


# xy = Array{Float64,2}
# xyz = Array{Float64,3}
# xyz = zeros()
# xyz[:,:,1] = xy

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

# Zf = individual_char(P,produc_t)
# A[1] = Matrix{Float64}
# typeof(A[:,1])  
# typeof(Z_test)
# A[:,:,1] = Z_test
# Z_test = individual_char()
# Z = individual_char()
# Z_fin = Matrix{Float64}
# Z_again = hcat(Z, Z_test)
# Z_fin = zeros(31,2)
# xy =[Z Z_test]
# Z_fin[:,1] = Z_test
# A = [f(a[i], b[j]) for i = 1:2, j = 1:3]
# 2×3 Matrix{Int64}:
# produc_t = precompute(carch_ch_mt,years,produc_t)
# prices = mz[Int64.(produc_t[1])]
# mz = Matrix()
# Array{Matrix{Float64},1}()
# produc_t = Array{Array{Float64,1},1}()
# xyz1 = zeros(Array{Float64,1},(31,1))
# produc_t[1]

# a = [1,2]
# push!(produc_t,[1.0,2.0,3.0])

# function make_price()
#     price_vec = zeros(6103,31)
#     for  i = 1:31
#         for j = 1:1
#             z[i, j] = mz[aProductID[i],j].*ind_characteristics'
#             price_vec[i,j] = car_ch[(car_ch[:,"Model_id"].==i).&(car_ch[:,"Year"].==j),:][:,"price"]
#         end

#     end
#     return price_vec
# end

# simdist'



# price_vector = make_price()

# mz = car_ch[:,"price"]

function value(lambda,t,Z)   ##Evaluating the value of individual characteristics,and prices u
    mu = Array{Float64,2}
    l = 1 ##Number of nonlinear attributes
    mu = lambda.*Z[t]
    mu_exp = exp.(mu)
    return mu_exp
end

# mu = value(P.lambda,1,Zf)
# Int64.(produc_t[1])
# R.delta_iia[327]
# R.delta_iia[Int64.(produc_t[1])]
# ev = exp.(R.delta_iia[Int64.(produc_t[1])])#.*mu
# ev./(1 .+ ev)
# test = exp.(delta_iia[Int64.(produc_t[1])]).*value(0.6,1)
# ms = test./(1 .+ sum_cols)
# sum_cols = sum(test,dims=1)
# mean(ms, dims=2)
# mat2 = ms*ms'./100
# ms.*(1 .- ms)
# Diagonal(mean(ms.*(1 .- ms),))

function demand(mu,t,delta,produc_t) ##Creating the Demand Function
    ev = exp.(delta[Int64.(produc_t[t])]).*mu
    market_share = ev./(1 .+ ev)
    #println("market share is ", market_share)
    shat = mean(market_share,dims=2)
    return market_share, shat

end

function jacobian(ms)

    square = ms*(1 .- ms)'
    mat2 = ms*ms'./100
    jacob = square./100 *I .-  (mat2 - mat2*I)/100
    return jacob
end





function inverse(P::Param,R::Results,eps0,eps1,share_ch,T_thisrun,shat,delta_iia,mu,produc_t,Z)
    for t = 1:T_thisrun
        value(P.lambda,t,Z)
        rowid = produc_t[t]
        f = 100
        res_norm::Array{Float64,1} = [0]
        i = 1
        jacob = Matrix{Float64}
        while norm(f) > eps1
            #jacob = 0
            ms, shat= demand(mu,T_thisrun,delta_iia,produc_t)
            println(i)
            #println("shat is",shat)
            jacob = jacobian(ms)
            # println(P.share_ch[Int64.(rowid)])
            # println(log.(shat))
            # println(log.(P.share_ch[Int64.(rowid)]))
            f = log.(P.share_ch[Int64.(rowid)]) .- log.(shat)
            # print(f)
            R.delta_iia[Int64.(rowid)] = R.delta_iia[Int64.(rowid)] .+ f
            println(norm(f))
            append!(res_norm,norm(f))
            print("res norm really is ", res_norm)
            i += 1
        end 

        while norm(f) > eps0
            #jacobian = 1
            ms,shat = demand(mu,T_thisrun,delta_iia,produc_t)
            println("shat2 is ", shat)
            f = log.(P.share_ch[Int64.(rowid)])-log.(shat)
            jacob = jacobian(ms)
            println(" Why is this Na ", R.delta_iia[Int64.(rowid)])
            println(" Why is this f Na ", f)
            R.delta_iia[Int64.(rowid)] = R.delta_iia[Int64.(rowid)] .+ inv(jacob./shat)*f
            #println(norm(f))
            append!(res_norm,norm(f))  
            println("res norm really really is ", res_norm)                         
        end     
        return res_norm   
    end
    println("res norm is ", res_norm)
    return res_norm
end 










function routine_overall(P::Param,R::Results,t)
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
    println("shat is ", shat)
    jacob = jacobian(market_share)
    #print(jacob)
    res_norm_init = inverse(P,R,P.eps0,P.eps1, P.share_ch,t,shat,R.delta_iia,mu,produc_t,Zf)
   #print(res_norm_init)
   println("res norm init ", res_norm_init)
    return res_norm_init
end





