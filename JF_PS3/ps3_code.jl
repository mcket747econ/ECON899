using Parameters, DataFrames, StatFiles, Random, Distributions, Optim, LinearAlgebra

cd("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3")

ind_characteristics = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Simulated_type_distribution.dta"))
car_characteristics = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Car_demand_characteristics_spec1.dta"))
iv_specification = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Car_demand_iv_spec1.dta"))


@with_kw struct parameters
    lambda_p::Float64 = 0.6
    R::Int64 = 79









end





function overall_demand_j()
    dem = exp(alpha*car_characteristics[,"price"] + beta*car_characteristics[:,"not_price"] + lambda*ind_characteristics[i]*car_chracteristics[:,"price"] )/(1+ sum_j(exp(alpha*car_characteristics[:,"price"] + beta*car_characteristics[:,"not_price"]+ lambda*ind_characteristics[i]*car_characteristics[:,"price"])))
    








end 
