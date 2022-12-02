using Parameters, DataFrames, StatFiles, Random, Distributions, Optim, LinearAlgebra

cd("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3")

ind_characteristics = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Simulated_type_distribution.dta"))
car_ch = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Car_demand_characteristics_spec1.dta"))
iv_specification = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Car_demand_iv_spec1.dta"))






@with_kw struct Param
    T::Int64 = length(unique(car_ch[:,"Year"]))  ##The number of markets(Years)
    delta_iia::Array{Float64,1} = Vector(car_ch[:,"delta_iia"]) ##the iia delta values
    simdist::Array{Float64,1} = Vector(ind_characteristics[:,"Var1"]) ##Vector containing the simulated individual characteristics
    mz::Array{Float64,1} = Vector(car_ch[:,"price"]) ##Vector containing prices for all models
    lambda::Float64 = 0.6
    eps1::Float64 = 1 
    eps0::Float64 = 10e-12
    years::Array{Float64,1}= sort(unique(car_ch[:,"Year"]))








end

@with_kw mutable struct Results






end 


sort(years)
@with_kw function Initialize()
    P = Param()


    return P 
end 



price_vector = zeros(6103,31 )
z = zeros(31,1)
share_ch = Vector(car_ch[:,"share"])

carch_ch_mt = Matrix(car_ch)
mz = Vector(car_ch[:,"price"])
mz1 = Vector(car_ch[:,"price"])

mz_test = Array{Float64,2}
mz_test = zeros(6103,2)
mz_test[:,1] = mz1
mz_test[:,2] = x
prices = car_p[xyz1]

Random.seed!(1234)
x = rand(Uniform(10,20),6103)
# xyz = getindex(carch_ch_mt[:,2],1997)
# xyz1 = findall(x->x==1985,carch_ch_mt[:,2])
# xyz2 = findall(x->x==1985,carch_ch_mt[:,2])
# xyz3 = findall(x->x==1985,carch_ch_mt[:,2])

xy = [xyz1,xyz2,xyz3]
xy[1]
years = unique(car_ch[:,"Year"])

# Creating a product_id_code
function precompute(carch_ch_mt,years,produc_t)
    produc_1 = 0
    for i = 1:31
        produc_1 = findall(x->x==years[i],carch_ch_mt[:,2])
        push!(produc_t,produc_1)

    end
    return produc_t
end



# xy = Array{Float64,2}
# xyz = Array{Float64,3}
# xyz = zeros()
# xyz[:,:,1] = xy

function individual_char()
    Zf = Array{Array{Float64,2},1}() ##Structure that will contain all the prices associated with products in a particular market(Year)
    Z = zeros(31,1)
    for i=1:31
        for j=1:1
            Z = mz[Int64.(produc_t[i])]*simdist' 
            push!(Zf,Z)
            #Zf[i,j] = Z

        end
    end
    return Zf
end

A = zeros(31)
A[1] = Matrix{Float64}
typeof(A[:,1])  
typeof(Z_test)
A[:,:,1] = Z_test
Z_test = individual_char()
Z = individual_char()
Z_fin = Matrix{Float64}
Z_again = hcat(Z, Z_test)
Z_fin = zeros(31,2)
xy =[Z Z_test]
Z_fin[:,1] = Z_test
A = [f(a[i], b[j]) for i = 1:2, j = 1:3]
2×3 Matrix{Int64}:
produc_t = precompute(carch_ch_mt,years,produc_t)
prices = mz[Int64.(produc_t[1])]
mz = Matrix()
Array{Matrix{Float64},1}()
produc_t = Array{Array{Float64,1},1}()
produc_t[1]

a = [1,2]
push!(produc_t,[1.0,2.0,3.0])

function make_price()
    price_vec = zeros(6103,31)
    for  i = 1:31
        for j = 1:1
            z[i, j] = mz[aProductID[i],j].*ind_characteristics'
            price_vec[i,j] = car_ch[(car_ch[:,"Model_id"].==i).&(car_ch[:,"Year"].==j),:][:,"price"]
        end

    end
    return price_vec
end

simdist'



price_vector = make_price()

mz = car_ch[:,"price"]

function value(lambda,t)   ##Evaluating the value of individual characteristics,and prices u
    mu = Array{Float64,2}
    l = 1 ##Number of nonlinear attributes
    mu = lambda.*Z[t]
    mu_exp = exp.(mu)
    return mu_exp
end


test = exp.(delta_iia[Int64.(produc_t[1])]).*value(0.6,1)
ms = test./(1 .+ sum_cols)
sum_cols = sum(test,dims=1)
mean(ms, dims=2)
mat2 = ms*ms'./100
ms.*(1 .- ms)
Diagonal(mean(ms.*(1 .- ms),))

function demand(mu,t,jacobian,delta,produc_t) ##Creating the Demand Function
    ev = exp.(delta_iia[Int64.(produc_t[t])]).*mu
    market_share = ev./(1.+ ev)
    shat = mean(market_share,dims=2)
    return market_share, shat

end
function jacobian(ms,)

    square = ms*(1 .- ms)'
    jacob = square./100 *I .-  (mat2 - mat2*I)/100
    return jacob

end



function inverse(epsilon::Float64,share_ch)
    for t = 1:T
        value(lambda,t)
        rowid = produc_t[t]
        f = 100
        res_norm = zeros(T,500)
        i = 1
        if norm(f) > eps1
            jacobian = 0
            ms, shat= demand(mu,t,jacobian,delta,produc_t)
            jacob = jacobian(ms,)
            f = log(share_ch[rowid])-log(shat)
            delta_iia[rowid] = delta_iia[rowid] + f
            res_norm[t, i] = norm(f)
            i += 1
        else 
            if norm(f) > eps0
                jacobian = 1
                ms,shat = demand(mu,t,jacobian,delta_produc_t)
                f = log(share_ch[rowid])-log(shat)
                jacobian = jacobian(ms,)
                delta_iia[rowid] .= delta_iia[rowid] .+ invert(jacobian./shat).*f
                res_norm[t] = norm(f)                                   
            end

        end
    end
    return res_norm, 
end 





 


square = ms*(1 .- ms)'

(mat2 - mat2*I)
1 .- I
mat2*I
mean(square,dims=2)
exp.(delta_iia[Int64.(produc_t[1])])
"]
car_ch[:,Model][5]

@with_kw struct parameters
    lambda_p::Float64 = 0.6
    R::Int64 = 79
    length_carch::Int64 = length(car_ch)
    delta_0::Array{Float64,1} = 
    










end

@with_kw mutable struct results
    demand::Array{Float64,2} 
    








end

function initialize()





end





function overall_demand_j(delta_0)
    dem = exp(alpha*car_characteristics[,"price"] + beta*car_characteristics[:,"not_price"] + lambda*ind_characteristics[i]*car_chracteristics[:,"price"] )/(1+ sum_j(exp(alpha*car_characteristics[:,"price"]
     + beta*car_characteristics[:,"not_price"]+ lambda*ind_characteristics[i]*car_characteristics[:,"price"])))
    







    return dem 
end 



function inversion(P)

    delta_0 = 
     for i = 1:NJ
        for i = 1:NT
            delta_test = P.delta_0[i,j] + (log(car_characteristics[i,"share"])-log(mean(demand(P.delta_0[i,j])))) 












     end





end


function pre_compute(P)
    
    Z_vec = zeros(P.T, P.Nb)
    for i = 0:P.T
        Z_vec = 
        for j=0:
    
    Z_vec[i,j] = Car_demand_

        




end


function routine_overall()






end
