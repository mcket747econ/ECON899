using Parameters, DataFrames, StatFiles, Random, Distributions, Optim, LinearAlgebra

cd("C:/Users/Rafeh/Documents/GitHub/ECON899/PS5_JF")

@with_kw struct Param
    δ::Float64 =0.1
    β::Float64 = 1/1.05
    α::Float64 = 0.06
    a::Float64 = 40.0
    b::Float64 = 10.0
    q_min::Int64 = 0
    q_max::Int64 = 45
    q_step::Int64 = 5
    q_grid::Array{Float64,1} = collect(q_min:q_step:q_max)
    nq::Int64 = length(q_grid)
    T::Int64 = 25
    xbar::Float64 = 20.0
    θ::Float64 = 0.5
end

# Store results when iterating 
# depreciation rate, capacity constraints
mutable struct Results
    δ::Float64
    α::Float64

    ## Investments 
    xopt_1::Array{Float64}
    xopt_2::Array{Float64}

    ## Value Functions 
    vf_1::Array{Float64}
    vf_2::Array{Float64}

    ## Profits
    πopt_1::Array{Float64}
    πopt_2::Array{Float64}

    ## Quantities
    qopt_1::Array{Float64}
    qopt_2::Array{Float64}
    ## Prices
    popt_1::Array{Float64}
    popt_2::Array{Float64}
    μ::Array{Float64}
end



function Initialize()
    P = Param()

    α = P.θ/((1-P.δ-P.θ)*P.xbar)
    ## Investments 
    xopt_1::Array{Float64} = zeros(P.nq, P.nq)
    xopt_2::Array{Float64} = zeros(P.nq, P.nq)

    ## Value Functions 
    vf_1::Array{Float64} = zeros(P.nq, P.nq)
    vf_2::Array{Float64} = zeros(P.nq, P.nq)

    ## Profits
    πopt_1::Array{Float64} = zeros(P.nq, P.nq)
    πopt_2::Array{Float64} = zeros(P.nq, P.nq) 

    ## Quantities
    qopt_1::Array{Float64} = zeros(P.nq, P.nq)
    qopt_2::Array{Float64} = zeros(P.nq, P.nq)
    ## Quantities
    popt_1::Array{Float64} = zeros(P.nq, P.nq)
    popt_2::Array{Float64} = zeros(P.nq, P.nq)
    μ = ones(P.nq, P.nq) 
    μ = μ./ (P.nq*P.nq)
    R = Results(P.δ,α,xopt_1,xopt_2,vf_1,vf_2,πopt_1,πopt_2,qopt_1,qopt_2,popt_1,popt_2,μ)
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

function probs(P::Param,vf,q,q_n,x1)
    @unpack q_min,q_max,q_step = P
    @unpack α,δ = R
    pr_inc = (1-δ)*α*x1 / (1+α*x1)
    pr_same = (1-δ)/(1+α*x1) + δ*α*x1 /(1+α*x1)
    pr_dec = (δ)/(1+α*x1)
    Δq = (q_n-q)/q_step
    if q<=q_min
        return [0.0,pr_same+pr_dec,pr_inc][max(1,min(3,trunc(Int64,Δq+2)))]
    elseif q>=q_max
        return [pr_dec,pr_same+pr_inc,0.0][max(1,min(3,trunc(Int64,Δq+2)))]
    else
        return [pr_dec,pr_same,pr_inc][max(1,min(3,trunc(Int64,Δq+2)))]
    end
    return 0.0
end
function cournot(P::Param,q_1,q_2,tol::Float64=1e-6)
    @unpack a = P
    qopt_1,qopt_1n = q_1,q_1
    qopt_2,qopt_2n = q_2,q_2
    diff = Inf 
    while (diff>tol)
        ## Bisection
        qopt_1n = max(0,min(q_1,0.5*(P.a-qopt_2)))
        qopt_2n = max(0,min(q_2,0.5*(P.a-qopt_1)))
        diff = max(norm(qopt_1n - qopt_1), norm(qopt_2n - qopt_2))
        qopt_1 = qopt_1n
        qopt_2 = qopt_2n
    end
    return qopt_1, qopt_2
end

function static_results(R::Results, P::Param)
    @unpack nq,q_grid = P
    @unpack qopt_2, qopt_1,popt_1,popt_2,πopt_1,πopt_2 = R
    for i = 1:nq
        for  j = 1:nq
            qopt_1[i,j], qopt_2[i,j] = cournot(P,q_grid[i], q_grid[j])
            pr = price(P,qopt_1[i,j] + qopt_2[i,j])
            popt_1[i,j] = pr
            popt_2[i,j] = pr
            πopt_1[i,j] = popt_1[i,j] * qopt_1[i,j]
            πopt_2[i,j] = popt_2[i,j] * qopt_2[i,j]
        end
    end
    return R
end


# written from the perspective of firm 1
function W(P::Param,R::Results,vf,q_1::Float64, q_2::Float64, x_2::Float64,δ, α)
    # get current index of firm 1
    @unpack q_min, q_step, nq,q_grid = P
    i = min(max(Int((q_1 -q_min)/q_step + 1),1),10)
    res = 0.0
    for j = 1:nq
        Δq = q_grid[j]-q_2
        # println(i,j)
        res += vf[i,j]*probs(P,R,q_2,q_grid[j], x_2)
    end
    return res
end


function bellman( P::Param,R::Results)
    @unpack δ, α, xopt_1, xopt_2, πopt_1,πopt_2, vf_1, vf_2 = R
    @unpack β, nq, q_grid,q_step = P

    # initialize updated guess
    x_1_n,x_2_n,vf_1_n,vf_2_n = zeros(nq, nq), zeros(nq, nq), zeros(nq, nq), zeros(nq, nq)
    for i = 1:nq
        for  j = 1:nq
            # quantities in grid
            q_1,q_2 = q_grid[i], q_grid[j]
            ## Firm 2's investment 
            # println(size(xopt_2))
            # println(i,j)
            x_2 = xopt_2[i, j]
            ## W under increasing/decreasing 
            W_const = W(P,R,R.vf_1,q_1,q_2, x_2,δ,α)
            W_dec = min(W_const,W(P,R,R.vf_1,q_1-q_step,q_2, x_2,δ,α))
            W_inc = max(W_const,W(P,R,R.vf_1,q_1+q_step,q_2, x_2,δ,α))

            ## Get best response from FOC
            if i == 1
                best_response = (-1 + sqrt(β*α*((1-δ)*(W_inc - W_const)) ))/α
            elseif i == nq
                best_response = (-1 + sqrt(β*α*(δ*(W_const - W_dec))))/α
            else
                best_response = (-1 + sqrt(β*α*((1-δ)*(W_inc - W_const) + δ*(W_const - W_dec)) ))/α
            end 
            ## Next period investment at i,j 
            x_1_n[i,j] = max(0, best_response)

            ## Subtract investment from π
            vf_1_n[i,j] = πopt_1[i,j] - x_1_n[i,j]
            for i_p = 1:nq
                vf_1_n[i,j] += β*W(P,R,R.vf_1,q_grid[i_p],q_2, x_2,δ,α)*probs(P,R,q_1, q_grid[i_p], x_1_n[i,j])
            end

            # solve firm 2 problem
            vf2_tmp = vf_2'
            x_1 = xopt_1[i, j]
            W_const = W(P,R,R.vf_1,q_2,q_1, x_1,δ,α)
            W_dec = min(W(P,R,vf2_tmp,q_2-q_step,q_1, x_1,δ, α), W_const)
            W_inc = max(W(P,R,vf2_tmp,q_2+q_step,q_1, x_1,δ, α), W_const)
            ## Get best response from FOC
            if i == 1
                best_response = (-1 + sqrt(β*α*((1-δ)*(W_inc - W_const)) ))/α
            elseif i == nq
                best_response = (-1 + sqrt(β*α*(δ*(W_const - W_dec))))/α
            else
                best_response = (-1 + sqrt(β*α*((1-δ)*(W_inc - W_const) + δ*(W_const - W_dec)) ))/α
            end 
            x_2_n[i,j] = max(0, best_response)
            vf_2_n[i,j] = πopt_2[i,j] - x_2_n[i,j]
            for j_p = 1:nq
                vf_2_n[i,j] += β * W(P,R,R.vf_1,q_grid[j_p],q_1,x_1,δ,α)*probs(P,R,q_2, q_grid[j_p], x_2_n[i,j])
            end     
        end
    end
    x_1_n, x_2_n, vf_1_n, vf_2_n
end

function solve_firm_prob(R::Results, P::Param)
    @unpack xopt_1, xopt_2,vf_1,vf_2 = R
    i, maxiter, err = 1, 10000, 100
    i = 1
    
    while (err > 1e-4) & (maxiter > i)
        x_1_n, x_2_n, vf_1_n, vf_2_n = bellman(P,R)
        err = max(norm(x_1_n - xopt_1), norm(x_2_n - xopt_2), norm(vf_1_n - vf_1), norm(vf_2_n - vf_2))
        xopt_1 = x_1_n
        xopt_2 = x_2_n
        vf_1 = vf_1_n
        vf_2 = vf_2_n

        i+=1
    end
    println("Converged in ", i)
end


function distribution(R::Results, P::Param)
    @unpack nq, q_grid = Param()
    μ_next = zeros(nq, nq)
    for i = 1:nq, j = 1:nq
        for i_p = 1:nq, j_p = 1:nq
            pr_1 = probs(P,R.vf_1,q_grid[i], q_grid[i_p], R.xopt_1[i,j])
            pr_2 = probs(P,R.vf_1,q_grid[j], q_grid[j_p], R.xopt_1[i,j])
            μ_next[i_p, j_p] += R.μ[i,j] * pr_1 * pr_2
        end
    end

    return μ_next
end

function solve_distrib(R::Results, P::Param; verbose::Bool = false)
    i, maxiter, err = 1, 10000, 100
    while (err > 1e-4) & (maxiter > i)
        μ_next = distribution(R, P)
        err = norm(R.μ - μ_next)
        R.μ = copy(μ_next)
        i+=1
    end

    println("Done, converged in ", i, " iterations.")
end

function Solve_model(P::Param,R::Results,δ::Float64)
    print("Computing quantities")
    R = static_results(R, P)

    print("Solving firm policy")
    solve_firm_prob(R, P)

    print("Stationary distribution")
    solve_distrib(R, P)
    return R
end



Ropt = Solve_model(P,R,0.1)
## Plotting


using Plots

p1 = plot(1:P.nq, 1:P.nq, Ropt.qopt_1, st=:surface, legend = false);
title!("q_star_1");
savefig("q_static_1.png")

p2 = plot(1:P.nq, 1:P.nq, Ropt.qopt_2, st=:surface, legend = false);
title!("q_star_2");
savefig("q_static_2.png")

p5 = plot(1:P.nq, 1:P.nq, Ropt.πopt_1, st=:surface, legend = false);
title!("π_1");
savefig("Profit_1.png")

p6 = plot(1:P.nq, 1:P.nq, Ropt.πopt_2, st=:surface, legend = false);
title!("π_2");
savefig("Profit_2.png")

p1 = plot(1:P.nq, 1:P.nq, Ropt.xopt_1, st=:surface, legend = false);
title!("Policy Function (1)");
savefig("xopt_1.png")

p2 = plot(1:P.nq, 1:P.nq, Ropt.xopt_2, st=:surface, legend = false);
title!("Policy Function (2)");
savefig("xopt_2.png")

p3 = plot(1:P.nq, 1:P.nq, Ropt.vf_1, st=:surface, legend = false);
title!("VF 1");
savefig("VF_1.png")

p4 = plot(1:P.nq, 1:P.nq, Ropt.vf_2, st=:surface, legend = false);
title!("VF 2");
savefig("VF_2.png")

p5 = plot(1:P.nq, 1:P.nq, Ropt.μ, st=:surface, legend = false);
title!("Stationary Distribution");
savefig("distribution.png")

