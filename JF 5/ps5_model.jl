using Parameters, DataFrames, StatFiles, Random, Distributions, Optim, LinearAlgebra

cd("C:/Users/Rafeh/Documents/GitHub/ECON899/PS5_JF")

@with_kw struct Param
    δ =0.1
    β = 1/1.05
    α = 0.06
    a = 40
    b = 10
    q_grid::Array{Float64,1} = 0:5:45
    nq::Int64 = length(q_grid)
    T::Int64 = 25
    q_step::Int64 = 5
end

# Store results when iterating 
# depreciation rate, capacity constraints
mutable struct Results
    δ::Float64
    α::Float64

    ## Investments 
    x_1::Array{Float64}
    x_2::Array{Float64}

    ## Value Functions 
    vf_1::Array{Float64}
    vf_2::Array{Float64}

    ## Profits
    π_1::Array{Float64}
    π_2::Array{Float64}

    ## Quantities
    q_1::Array{Float64}
    q_2::Array{Float64}
    Results() = new()
end



function Initialize()
    P = Param()
    R = Results()
    return P, R
end 

P,R = Initialize()

Random.seed!(1234)
function demand(P::Param,price::Float64)
    return P.a-P.b*price
end
## Q = a - bP => P = (a - Q)/b

function price(P::Param,quantity::Float64)
    return (P.a-quantity)/P.b
end

function probs(Δq,x1)
    pr = 0.0
    if Δq = 1
        pr = (1-δ)*α*x1 / (1+α*x1)
    elseif Δq == 0
        pr = (1-δ)/(1+α*x1) + δ*α*x1 /(1+α*x1)
    elseif Δq==-1
        pr = (δ)/(1+α*x1)
    end
    return pr
end
function cournot(P::Param,q_1,q_2)
    @unpack a = P
    qopt_1,qopt_1n = q_1,q_1
    qopt_2,qopt_2n = q_2,q_2
    diff = Inf 
    while (diff>tol)
        ## Bisection
        qopt_1n = max(0,min(q_1,0.5*(P.a-qopt_2)))
        qopt_2n = max(0,min(q_2,0.5*(P.a-qopt_1)))
        diff = max(norm(qopt_1n - qopt), norm(qopt_2n - qopt_2))

        qopt_1 = qopt_1n
        qopt_2 = qopt_2n
    end
    return qopt_1, qopt_2
end

function static_results(R::Results, P::Primitives)
    @unpack nq = P
    @unpack q_2, q_1 = R
    for i = 1:nq
        for  j = 1:nq
            q_1[i,j], q_2[i,j] = cournot(q_grid[i], q_grid[j], P)
            p_star_1[i,j] = inv_demand(q_1[i,j] + q_2[i,j], P)
            p_star_2[i,j] = inv_demand(q_1[i,j] + q_2[i,j], P)
            π_1[i,j] = p_1[i,j] * q_1[i,j]
            π_2[i,j] = p_2[i,j] * q_2[i,j]
        end
    end
end


# written from the perspective of firm 2
function W(q_1::Float64, q_2::Float64, x_2::Float64, δ::Float64, α::Float64, vf::Array{Float64}, P::Param)
    # get current index of firm 1
    @unpack q_min, q_step, nq = P
    i = Int((q_1 -q_min)/q_step + 1)
    res = 0.0
    for j = 1:nq
        Δq = q_grid[j]-q_2
        res += vf[i,j]*pr(Δq, x_2, δ, α, P)
    end

    return res
end


function solve_firm_problem!(R::Results, P::Param, γ::Float64 = 1.0)
    @unpack x_1 = R
    while (err > 1e-4) & (maxiter > i)
        x_1_n, x_2_n, vf_1_next, vf_2_next = bellman(R, P)
        err = max(norm(x_1_n - R.x_pf_1), norm(x_2_n - R.x_pf_2), norm(vf_1_next - R.vf_1), norm(vf_2_next - R.vf_2))
        x_pf_1 = (1-γ)*x_pf_1 + γ*copy(x_1_n)
        x_pf_2 = (1-γ)*x_pf_2 + γ*copy(x_2_n)
        vf_1 = (1-γ)*vf_1 + γ*copy(vf_1_next)
        vf_2 = (1-γ)*vf_2 + γ*copy(vf_2_next)

        i+=1
    end

    println("Done, converged in ", i, " iterations.")
end

## FOC is $
function bellman(q_1,q_2,tol=1e-8)
    π = profit(q_1,q_2)
    sum = 0

    for x1 in x1_grid
        if x1>0
            sum+=-x1
            for q1 in q_grid
                for q2 in q_grid
                probs()
                end
            end
        end
    end
    return produc_t
end

