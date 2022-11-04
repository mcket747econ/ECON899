using Random, Distributions, DataFrames, Parameters, Optim, Plots

@with_kw struct parameters
    rho_0::Float64 = 0.5
    sigma_0::Float64 = 1
    x_0::Float64 = 0
    T::Int64 = 200
    l::Int64 = 2
    H::Int64 = 10






end

@with_kw mutable struct results
    x_t::Array{Float64,2}
    e::Array{Float64,2}
    M::Array{Float64,1}
    # shocks::Array{Float64,1}






end

function initialize()
    P = parameters()
    x_t = zeros(P.T,P.H)
    e = rand(Normal(0,(P.sigma_0)^2),P.T,P.H)
    # shocks = zeros(P.T)
    M = zeros(2)
    R = results(x_t,e,M)









    return P, R
end

Random.seed!(1234)
P,R = initialize()

function simulate(P,e)
    #set.seed(123)
    x_t = zeros(P.T,P.H)
    x_t[1,:] = P.rho_0*P.x_0 .+ e[1,:]
    for i in 2:P.T
        x_t[i,:] = P.rho_0*x_t[i-1,:] .+ e[i,:]
    end

    return x_t

end

###Initial_simulation
R.x_t = simulate(P,R.e)
R.M = [mean(R.x_t), sum((R.x_t .- mean(R.x_t)).*((R.x_t .- mean(R.x_t))/P.T))]
