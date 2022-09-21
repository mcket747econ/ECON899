# using Pkg; Pkg.add(["Plots","Parameters"])
using Plots, Parameters
@with_kw struct Primitives
    β::Float64 = 0.9932 #Discount Rate
    b::Float64 = 0.5 #Unemployment Benefits
    α::Float64 = 1.5
    price::Array{Float64,1} = collect(range(0,length=1000,stop=1))
    S::Array{Float64,1} = [1, 0.5]
    Π::Array{Float64,2} = [0.97 .03; 0.5 0.5]
    A:: Array{Float64,1}= collect(range(-2,length=1000,stop=5))
    ne::Int64 = length(S)
    na::Int64 = length(A)


end

mutable struct Results
    val_func:: Array{Float64,2}
    pol_func:: Array{Float64,2}
    mu:: Array{Float64,2}
    q:: Float64
end0


function Initialize()
    prim = Primitives()
    val_func = zeros(prim.na,prim.ne)
    pol_func = zeros(prim.na,prim.ne)
    mu = ones(prim.na, prim.ne)/(prim.na*prim.ne)## The division is in order to create a uniform distribution
    q = (1 + prim.β)/2
    res = Results(val_func, pol_func, mu, q)
    # V_iterate(prim,res)
    prim, res



end
function V_iterate(prim::Primitives, res::Results, tol::Float64 = 1e-3)
    error = 100
    n = 0
    while error>tol
        n+=1
        v_next = Bellman(prim,res)
        error = maximum(abs.((v_next-res.val_func))) #I assume this is the maximum of the array
        res.val_func = v_next
        println("Iteration number: ", n)
        println("Absolute difference: ", error)
    end
    print("Value Function Converged in ", n, " Iterations.")

end

function Bellman(prim::Primitives,res::Results)
    @unpack  Π, β, na, ne, price, α, b, A, S = prim
    v_next = zeros(na, ne)

    for s_index=1:ne
        ap_index_start = 1
        for a_index = 1:na
            v_0 = -1e10
            budget = S[s_index] + A[a_index]
            for ap_index=ap_index_start:na
                c = budget - A[ap_index]*res.q
                if c > 0
                    v = ((c^(1-α))-1)/(1-α) + β*sum(res.val_func[ap_index,:].*Π[s_index,:])
                    if v > v_0
                        v_0 = v
                        res.pol_func[a_index,s_index] = A[ap_index]
                        ap_index_start = ap_index
                    end
                end
            end
            v_next[a_index,s_index] = v_0
        end
    end
    v_next
end

function distribution(prim::Primitives, res::Results, tol::Float64 =-1e3)
    @unpack Π, β, na, ne, price, α, b, A, S = prim
    mu_next = zeros(na, ne)
    # mu_0 = ((A[0]-(-2))/(7))
    for a_index=1:na,s_index=1:ne
        ap = res.pol_func[a_index,s_index]
        ap_index = argmin(abs.(ap .- A))
        #s_p = res.pol_func(a_index,:)
        for sp_index = 1:ne
            mu_next[ap_index,sp_index] += res.mu[a_index, s_index] * Π[s_index, sp_index]
        end

    end
    mu_next
end


function dist_iterate(prim::Primitives, res::Results, tol::Float64 = 1e-3)
    error = 100
    n = 0
    while error>tol
        n+=1
        mu_next = distribution(prim,res)
        error = maximum(abs.(mu_next-res.mu))#I assume this is the maximum of the array
        res.mu = mu_next
        println("Iteration number: ", n)
        println("Absolute difference: ", error)
    end
    print("Mu Distribution Converged in ", n, " Iterations.")

end

function price_update(prim::Primitives, res::Results, tol::Float64 = 1e-3)
    aggregate_a = sum(res.mu[:,1] .* prim.A ) + sum(res.mu[:,2] .* prim.A ) #if it is positive: oversaving, negative: overborrowing.
    if abs(aggregate_a) > tol
        if aggregate_a > 0
            res.q = res.q + (1 - res.q)/2 * abs(aggregate_a)
        elseif aggregate_a < 0
            res.q = res.q - (res.q-prim.β)/2 * abs(aggregate_a)
        end
    end
end

prim, res = Initialize()


V_iterate(prim,res)
dist_iterate(prim,res)
price_update(prim, res)

function overall_iterate(prim::Primitives, res::Results, tol::Float64 = 1e-4)
    balance = 100
    x = 0
    while abs(balance) > tol
        x+=1
        V_iterate(prim,res)
        dist_iterate(prim,res)
        balance = aggregate_a = sum(res.mu[:,1] .* prim.A ) + sum(res.mu[:,2] .* prim.A )
        price_update(prim,res)
        println("Iteration number: ", x)
        println("Absolute difference: ", balance)
    end


end

overall_iterate(prim,res)

Plots.plot(prim.A, res.mu, title="Asset Distribution", label = ["Employed" "Unemployed"])
sum(res.mu[:,])
sum(res.mu[:,1])
sum(res.mu[:,2])


end
