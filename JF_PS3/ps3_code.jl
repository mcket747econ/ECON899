using Parameters, DataFrames, StatFiles, Random, Distributions, Optim, LinearAlgebra

cd("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3")

ind_characteristics = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Simulated_type_distribution.dta"))
car_ch = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Car_demand_characteristics_spec1.dta"))
iv_specification = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Car_demand_iv_spec1.dta"))

car_ar = convert(Array, car_ch)
Ma
car_ch[,]
price_vector = zeros(6103,31 )
z = zeros(31,1)

Array(car_ch)





end





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





price_vector = make_price()

mz = car_ch[:,"price"]


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










